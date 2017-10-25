#include <cmath>
#include "ThermalConvection.hpp"
#include "../Mara.hpp"
#include "../BoundaryConditions.hpp"
#include "../BlockDecomposition.hpp"
#include "../CartesianMeshGeometry.hpp"
#include "../CellCenteredFieldCT.hpp"
#include "../Checkpoint.hpp"
#include "../ConservationLaws.hpp"
#include "../FieldOperator.hpp"
#include "../IntercellFluxSchemes.hpp"
#include "../MeshData.hpp"
#include "../MeshOperator.hpp"
#include "../Problems.hpp"
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
        Cow::Array& A,
        MeshLocation location,
        MeshBoundary boundary,
        int axis,
        int numGuard) const override
    {
        switch (axis)
        {
            case 0: return periodic.apply (A, location, boundary, axis, numGuard);
            case 1: return periodic.apply (A, location, boundary, axis, numGuard);
            case 2:

            if (boundary == MeshBoundary::left)
                return reflecting.apply (A, location, boundary, axis, numGuard);

            for (int i = 0; i < A.size(0); ++i)
            {
                for (int j = 0; j < A.size(1); ++j)
                {
                    for (int k = A.size(2) - numGuard; k < A.size(2); ++k)
                    {
                        auto I = Index {{i, j, k - numGuard, 0, 0}};
                        auto z = geometry->coordinateAtIndex(I)[2];
                        auto P = idf (0, 0, z);

                        for (int q = 0; q < 5; ++q)
                        {
                            A(i, j, k, q) = P[q];
                        }
                    }
                }
            }
        }
    }

    bool isAxisPeriodic (int axis) override
    {
        switch (axis)
        {
            case 0: return true;
            case 1: return true;
            case 2: return false;
            default: throw std::logic_error ("ThermalConvectionBoundaryCondition");
        }
    }

    ReflectingBoundaryCondition reflecting;
    PeriodicBoundaryCondition periodic;
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
    user["tfinal"]  = 16.0;
    user["serial"]  = false;
    user["cpi"]     = 0.25;
    user["cpf"]     = "single"; // or multiple
    user["tsi"]     = 0.1;
    user["cfl"]     = 0.5;
    user["plm"]     = 2.0;
    user["N"]       = 16;
    user["aspect"]  = 1;
    user["dims"]    = 2;
    user["gamma"]   = 5. / 3;


    // Gravitational source terms, heating, and initial data function
    // ------------------------------------------------------------------------
    const double q0 = 2.0; // heating hate
    const double d0 = 1.0; // density at base of atmosphere
    const double g0 = 1.0; // gravitational acceleration at base
    const double z0 = 1.0; // gravity fall-off distance
    const double h = 0.25; // pressure scale height
    const double cs2 = h * g0; // nominal, 'isothermal' sound speed at base

    auto sourceTermsFunction = [&] (double x, double y, double z, StateArray p)
    {
        const double t = status.simulationTime;
        const double P = 1; // time period

        const double g = g0 / (1.0 + z / z0);
        const double rho = p[0];
        const double vel = p[3];
        auto S = StateArray();
        S[3] = -rho * g;
        S[4] = -rho * g * vel;

        double qdot = q0 * std::exp (-(y * y + z * z) / 0.01);
        S[4] += qdot * (t < P);

        return S;
    };

    auto initialData = [&] (double x, double y, double z) -> std::vector<double>
    {
        const double rho = d0 * std::pow (1 + z / z0, -z0 / h);
        const double pre = cs2 * rho;
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
    auto fs = std::make_shared<MethodOfLinesPlm>();
    auto mg = std::shared_ptr<MeshGeometry> (new CartesianMeshGeometry);
    auto cl = std::make_shared<NewtonianHydro>();


    // Set up grid shape. In 1D it's z-only. In 2D it's the x-z plane.
    // ------------------------------------------------------------------------
    auto dims = int (user["dims"]);
    int Nvert = int (user["N"]) * int (user["aspect"]);
    int Nhori = int (user["N"]);
    auto cs = Shape {{ Nhori, Nhori, Nvert, 1, 1 }};
    auto bs = Shape {{ 2, 2, 2, 0, 0 }};
    if (dims <= 2) { cs[0] = 1; bs[0] = 0; }
    if (dims <= 1) { cs[1] = 1; bs[1] = 0; }

    mg->setCellsShape (cs);
    mg->setLowerUpper ({{-0.5, -0.5, 0.0}}, {{0.5, 0.5, 1.0 * int (user["aspect"])}});

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




