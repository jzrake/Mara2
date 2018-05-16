#include <cmath>
#include <cassert>
#include "ThermalConvection.hpp"
#include "../Mara.hpp"
#include "../BoundaryConditions.hpp"
#include "../BlockDecomposition.hpp"
#include "../CartesianMeshGeometry.hpp"
#include "../SphericalMeshGeometry.hpp" //include spherical mesh geometry
#include "../CellCenteredFieldCT.hpp"
#include "../Checkpoint.hpp"
#include "../ConservationLaws.hpp"
#include "../FieldOperator.hpp"
#include "../IntercellFluxSchemes.hpp"
#include "../MeshData.hpp"
#include "../MeshOperator.hpp"
#include "../Problems.hpp"
#include "../RiemannSolvers.hpp"
#include "../SolutionSchemes.hpp"
#include "../TaskScheduler.hpp"
#include "../TimeSeriesManager.hpp"
#include "MPI.hpp"
#include "HDF5.hpp"
#include "Timer.hpp"
#include "Logger.hpp"

using namespace Cow;




/**
Boundary conditions for a vertical atmosphere. The x and y directions are
periodic. The floor is reflecting and the ceiling is currently set to the
value of the initial condition.
*/
class ThermalConvectionBoundaryCondition : public BoundaryCondition
{
public:
    ThermalConvectionBoundaryCondition (InitialDataFunction idf) : idf (idf)
    {
        reflecting.setConservationLaw (std::make_shared<NewtonianHydro>());
    }

    void setMeshGeometry (std::shared_ptr<MeshGeometry> geometryToUse) override
    {
        geometry = geometryToUse;
    }

    void apply (
      Cow::Array& A, //boundary condition array
      MeshLocation location, //enum class {vert, edge, face, cell} in mara.hpp
      MeshBoundary boundary, //either left or right, enum class
      int axis, //0 ,1 or 2
      int numGuard) const override //number of guardzones; overrides virtual method in BoundaryConditions class in Mara.hpp
      {
          switch (axis) // this will be iterated in the applySimple method in BoundaryConditions.cpp, which is called in MeshData as applyBoundaryCondition
          {
              case 0:
              if (boundary == MeshBoundary::left) return reflecting.apply (A, location, boundary, axis, numGuard); //this apply belongs to derived class ReflectingBoundaryCondition
              if (boundary == MeshBoundary::right) return reflecting.apply (A, location, boundary, axis, numGuard);

              /*
              for (int i = A.size(0) - numGuard; i < A.size(0); ++i) //looping over radial guardzones
              {
                  for (int j = 0; j < A.size(1); ++j)
                  {
                      for (int k = 0; k < A.size(2); ++k)
                      {
                          auto I = Index {{i - numGuard, j, k, 0, 0}};
                          auto r = geometry->coordinateAtIndex(I)[0];
//                        auto q = geometry->coordinateAtIndex(I)[1]; //for querying coordinates in q and p
//                        auto p = geometry->coordinateAtIndex(I)[2];

                          auto P = idf (r, 0, 0); //change this to r,p,q in multidimensional

                          for (int q = 0; q < 5; ++q)
                          {
                              A(i, j, k, q) = P[q]; //q is hydrodynamic variable
                          }
                      }
                  }
              }
              */
//asserting to see if branch hit in 1D
//            case 1: assert(false); return reflecting.apply (A, location, boundary, axis, numGuard); //uses the BoundaryCondition::periodic apply method; overloaded method
//            case 2: assert(false); return periodic.apply (A, location, boundary, axis, numGuard); //applies periodic condition to non-r axes

          }
      }


    bool isAxisPeriodic (int axis) override
    {
        switch (axis)
        {
            case 0: return false;
            case 1: return false;
            case 2: return true;
            default: throw std::logic_error ("ThermalConvectionBoundaryCondition");
        }
    }

    ReflectingBoundaryCondition reflecting;
    PeriodicBoundaryCondition periodic;
    OutflowBoundaryCondition outflow;
    InitialDataFunction idf;
    std::shared_ptr<MeshGeometry> geometry;
};




// ============================================================================
int ThermalConvectionProgram::run (int argc, const char* argv[])
{
    auto status = SimulationStatus();
    auto user = Variant::NamedValues();
    user["outdir"]  = "data";
    user["restart"] = "";
    user["tfinal"]  = 3000; //1100
    user["serial"]  = false;
    user["cpi"]     = 0.25;
    user["cpf"]     = "single"; // or multiple
    user["tsi"]     = 0.1;
    user["cfl"]     = 0.5;
    user["plm"]     = 2.0;
    user["N"]       = 1000; //16 default
    user["aspect"]  = 1;
    user["dims"]    = 2;
    user["gamma"]   = 5. / 3;
    user["cartesian"] = true; //setting spherical or cartesian, default to cartesian
    user["g0"]      = 0.2; // g0 = G*M
    user["rho0"]    = 1.0; // inner density
    user["h"]       = 0.5; // pressure scale height, tweak to get strong/weak gravity
    user["K"]       = 1.0; // entropy

    // Gravitational source terms, heating, and initial data function
    // ------------------------------------------------------------------------

    const double g0 = user["g0"];
    const double rho0 = user["rho0"];
    const double h = user["h"];
    const double cs2 = h * g0;
    const double gamma = user["gamma"];
    const double K = user["K"];

    auto sourceTermsFunction = [&] (double r, double q, double p, StateArray P)
    {
        const double dg = P[0];
        const double vr = P[1];
        const double vq = P[2];
        const double vp = P[3];
        const double pg = P[4];

        auto S = StateArray();
        S[0] = 0.0;
        S[1] = (2 * pg + dg * (vq * vq + vp * vp)) / r;
        S[2] = pg * std::tan (M_PI_2 - q) / r + dg * (vp * vp * std::tan (M_PI_2 - q) - vr * vq) / r;
        S[3] = -dg * vp * (vr + vq * std::tan (M_PI_2 - q)) / r;
        S[4] = 0.0;

        const double g = g0 * std::pow(r,-2.0);
        const double rho = P[0];

        S[1] += -rho * g;
        S[4] += -rho * g * vr;

        return S;
    };

    auto initialData = [&] (double r, double q, double p) -> std::vector<double>
    {
        //*********************************************************************
        //hard-coded to assume r_0 = 1.0;
        /* const double alpha = (gamma - 1.0) / gamma;
        const double beta = g0 / K * alpha;
        const double v0 = std::pow(rho0,(gamma - 1.0));
        const double c1 = v0 + beta / (2.0* alpha - 1.0);
        const double rho = std::pow((c1 * std::pow(r,-2.0*alpha) - beta / (2.0* alpha - 1.0) * (1.0 / r) ),(1.0 / (gamma -1.0)));
        */
        const double rho = 1.0;
        const double pre = cs2 * rho / gamma;
        return std::vector<double> {rho, 0, 0, 0, pre};
    };


    auto clineUser = Variant::fromCommandLine (argc, argv);
    auto chkptUser = Variant::NamedValues();

    if (clineUser.find ("restart") != clineUser.end())
    {
        chkptUser = H5::File (clineUser["restart"], "r").getGroup ("user").readNamedValues();
    }
    Variant::update (user, chkptUser);
    Variant::update (user, clineUser);

    auto logger    = std::make_shared<Logger>();
    auto writer    = std::make_shared<CheckpointWriter>();
    auto tseries   = std::make_shared<TimeSeriesManager>();
    auto scheduler = std::make_shared<TaskScheduler>();
    double timestepSize = 0.0;

    auto bd = std::shared_ptr<BlockDecomposition>();
    auto bc = std::shared_ptr<BoundaryCondition> (new ThermalConvectionBoundaryCondition (initialData));
    auto rs = std::shared_ptr<RiemannSolver> (new HllcNewtonianHydroRiemannSolver);
    auto fs = std::shared_ptr<IntercellFluxScheme> (new MethodOfLinesPlm);
    auto cl = std::make_shared<NewtonianHydro>();
/*
Gridshape setlowerupper and .... needs to be changed here, scope problems here
*/
    auto bs = Shape();

    // Set up grid shape. In 1D it's z-only. In 2D it's the x-z plane.
    // ------------------------------------------------------------------------
/*    if (user["cartesian"])
    {
    auto mg = std::shared_ptr<MeshGeometry> (new CartesianMeshGeometry);
    auto dims = int (user["dims"]);
    int Nvert = int (user["N"]) * int (user["aspect"]); //gridpoints in vertical
    int Nhori = int (user["N"]); //gridpoints in horizontal
    auto cs = Shape {{ Nhori, Nhori, Nvert, 1, 1 }}; //square base gridpoints
    bs = Shape {{ 2, 2, 2, 0, 0 }};
    if (dims <= 2) { cs[0] = 1; bs[0] = 0; }
    if (dims <= 1) { cs[1] = 1; bs[1] = 0; }
    mg->setCellsShape (cs);
    mg->setLowerUpper({{-0.5, -0.5, 0.0}}, {{0.5, 0.5, 1.0 * int (user["aspect"])}});
    }
    else //spherical part here
    {
    */
    auto mg = std::shared_ptr<MeshGeometry> (new SphericalMeshGeometry);
    bs = Shape {{ 2, 0, 0, 0, 0 }};
    mg->setCellsShape ({{ user["N"], 1, 1 }});
    mg->setLowerUpper ({{ 1.0, M_PI*0.5-0.1, 0}}, {{10.0, M_PI*0.5+0.1, 0.1}});
  //  }

    if (! user["serial"])
    {
        logger->setLogToNullUnless (MpiCommunicator::world().isThisMaster());
        bd = std::make_shared<BlockDecomposition> (mg, *logger);
        mg = bd->decompose();
        bc = bd->createBoundaryCondition (bc);
    }
    tseries->setLogger (logger);
    scheduler->setLogger (logger);

    auto ss = std::make_shared<MethodOfLinesTVD>();
    auto mo = std::make_shared<MeshOperator>();
    auto fo = std::make_shared<FieldOperator>();
    auto md = std::make_shared<MeshData> (mg->cellsShape(), bs, cl->getNumConserved());

    cl->setGammaLawIndex (double (user["gamma"]));
    fs->setPlmTheta (double (user["plm"]));
    fs->setRiemannSolver (rs);
    fo->setConservationLaw (cl);
    mo->setMeshGeometry (mg);
    ss->setSourceTermsFunction (sourceTermsFunction);
    ss->setMeshOperator (mo);
    ss->setFieldOperator (fo);
    ss->setBoundaryCondition (bc);
    ss->setRungeKuttaOrder (2);
    ss->setDisableFieldCT (false);
    ss->setIntercellFluxScheme (fs);

    md->setVelocityIndex (cl->getIndexFor (ConservationLaw::VariableType::velocity));
    bc->setMeshGeometry (mg);
    bc->setConservationLaw (cl);

    auto L         = mo->linearCellDimension();
    auto V         = mo->measure (MeshLocation::cell);
    auto P         = md->getPrimitive();
    auto advance   = [&] (double dt) { return ss->advance (*md, dt); };
    auto condition = [&] () { return status.simulationTime < double (user["tfinal"]); };
    auto timestep  = [&] () { return timestepSize; };

    auto taskRecomputeDt = [&] (SimulationStatus, int rep)
    {
        double localDt = double (user["cfl"]) * fo->getCourantTimestep (P, L);

        timestepSize = MpiCommunicator::world().minimum (localDt);
    };

    auto taskCheckpoint = [&] (SimulationStatus, int rep)
    {
        writer->writeCheckpoint (rep, status, *cl, *md, *mg, *logger);
    };

//lambda for computing primitives averaged over cell, rep is repetition
    auto taskTimeSeries = [&] (SimulationStatus, int rep)
    {
        //bd defined, volumeAverageOverPatches returned, skips val forloop
        auto volumeAverageOverPatches = [&] (std::vector<double> vals)
        {
            if (bd)
                return bd->volumeAverageOverPatches (vals);

            for (auto& val : vals)
                val /= mg->meshVolume();

            return vals;
        };

        auto volumeIntegrated = fo->volumeIntegratedDiagnostics (P, V);
        auto volumeAveraged = volumeAverageOverPatches (volumeIntegrated);
        auto fieldNames = cl->getDiagnosticNames();
        auto entry = Variant::NamedValues();

        for (unsigned int n = 0; n < fieldNames.size(); ++n) //loop to set entry['primitive--name'] = volumeAveragedvalue, unrestricted in size
        {
            entry[fieldNames[n]] = volumeAveraged[n];
        }
        /*entry['total_entropy'] = volumeAveraged[]
        auto volumeIntegrated = fo->volumeIntegratedDiagnostics (P, V); //add specific entropy to be integrated
        auto volumeAveraged = volumeAverageOverPatches (volumeIntegrated);

        for (auto it = entry.begin(); it != entry.end(); it++){ //cout-ing name value pair--testing here
            std::cout << it->first << "and" << it->second << std::endl;
        }*/
        tseries->append (status, entry);
    };

//schedule time series data; tasktimeseries callback; compute energy per cell and log, scheduler passed into maraMainLoop which is updated
    scheduler->schedule (taskTimeSeries, TaskScheduler::Recurrence (user["tsi"]), "time_series");
    scheduler->schedule (taskCheckpoint, TaskScheduler::Recurrence (user["cpi"]), "checkpoint");
    scheduler->schedule (taskRecomputeDt, TaskScheduler::Recurrence (0.0, 0.0, 16), "compute_dt");

    writer->setMeshDecomposition (bd);
    writer->setTimeSeriesManager (tseries);
    writer->setTaskScheduler     (scheduler);
    writer->setOutputDirectory   (user["outdir"]);
    writer->setFormat            (user["cpf"]);
    writer->setUserParameters    (user);
    writer->setFilenamePrefix    ("chkpt");

    status = SimulationStatus();
    status.totalCellsInMesh = mg->totalCellsInMesh();

    if (! user["restart"].empty())
    {
        writer->readCheckpoint (user["restart"], status, *cl, *md, *logger);
        scheduler->skipNext ("checkpoint");
    }
    else
    {
        md->assignPrimitive (mo->generate (initialData, MeshLocation::cell));
    }
    md->applyBoundaryCondition (*bc);
    taskRecomputeDt (status, 0);

    logger->log() << std::endl << user << std::endl;
    return maraMainLoop (status, timestep, condition, advance, *scheduler, *logger);
}
