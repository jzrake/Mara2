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
#define cot(x) std::tan (M_PI_2 - x)
using namespace Cow;




/**
Boundary conditions for a vertical atmosphere. The x and y directions are
periodic. The floor is reflecting and the ceiling is currently set to the
value of the initial condition.
*/
class ThermalConvectionBoundaryCondition : public BoundaryCondition
{
public:
    ThermalConvectionBoundaryCondition()
    {
        reflecting.setConservationLaw (std::make_shared<NewtonianHydro>());
    }

    void apply (Cow::Array& A, MeshLocation location, MeshBoundary boundary, int axis, int numGuard) const override
    {
        switch (axis)
        {
            case 0:
            {
                if (boundary == MeshBoundary::left)
                {
                    return reflecting.apply (A, location, boundary, axis, numGuard);
                }
                if (boundary == MeshBoundary::right)
                {
                    return outflow.apply (A, location, boundary, axis, numGuard);
                }
            }
            case 1: return periodic.apply (A, location, boundary, axis, numGuard);
            default: throw std::logic_error ("ThermalConvectionBoundaryCondition: not yet configured for 3D");
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
};




// ============================================================================
int ThermalConvectionProgram::run (int argc, const char* argv[])
{
    auto status = SimulationStatus();
    auto user = Variant::NamedValues();
    user["outdir"]  = "data";
    user["restart"] = "";
    user["tfinal"]  = 100.0;
    user["serial"]  = false;
    user["cpi"]     = 0.25;
    user["cpf"]     = "single"; // or multiple
    user["tsi"]     = 0.1;
    user["cfl"]     = 0.5;
    user["plm"]     = 2.0;
    user["Nr"]      = 64;
    user["Nt"]      = 64;
    user["gamma"]   = 5. / 3;
    user["g0"]      = 0.2; // g0 = G*M
    user["rho0"]    = 1.0; // inner density
    user["h"]       = 0.5; // pressure scale height, tweak to get strong/weak gravity
    user["K"]       = 1.0; // entropy
    user["uniform"] = false;



    auto clineUser = Variant::fromCommandLine (argc, argv);
    auto chkptUser = Variant::NamedValues();

    if (clineUser.find ("restart") != clineUser.end())
    {
        chkptUser = H5::File (clineUser["restart"], "r").getGroup ("user").readNamedValues();
    }
    Variant::update (user, chkptUser);
    Variant::update (user, clineUser);



    const double g0    = user["g0"];
    const double rho0  = user["rho0"];
    const double h     = user["h"];
    const double gamma = user["gamma"];
    const double K     = user["K"];
    const double q0    = user["q0"];



    // Gravitational source terms, heating, and initial data function
    // ------------------------------------------------------------------------
    auto sourceTermsFunction = [g0, h, q0] (double r, double q, double p, StateArray P)
    {
        const double dg = P[0];
        const double vr = P[1];
        const double vq = P[2];
        const double vp = P[3];
        const double pg = P[4];
        auto S = StateArray();


        S[0] = 0.0;
        S[1] = (2 * pg + dg * (vq * vq + vp * vp)) / r;
        S[2] = (pg * cot(q) + dg * (vp * vp * cot(q) - vr * vq)) / r;
        S[3] = -dg * vp * (vr + vq * cot(q)) / r;
        S[4] = 0.0;


        const double g = g0 * std::pow (r, -2.0);
        S[1] += -dg * g;
        S[4] += -dg * g * vr;


        // exponential heating
        S[4] += q0 * std::exp (-1.0 * r / h);


        return S;
    };

    auto initialData = [&user, g0, rho0, gamma, K] (double r, double q, double p) -> std::vector<double>
    {
        using std::pow;

        double delta = 0.0; // 1e-1 * rand() / RAND_MAX;

        if (! user["uniform"])
        {
            auto r0 = 1.0;
            auto a = (gamma - 1.0) / gamma;
            auto d0 = rho0;
            auto c = pow(d0, gamma - 1) - g0 * a / K / r0;
            auto d = pow(g0 * a / K / r + c, 1. / (gamma - 1.));
            auto pre = K * pow(d, gamma);
            return {d * (1 + delta), 0, 0, 0, pre};
        }
        else
        {
            const double rho = 1.0;
            const double pre = 1.0;
            return std::vector<double> {rho, 0, 0, 0, pre};
        }
    };


    auto logger    = std::make_shared<Logger>();
    auto writer    = std::make_shared<CheckpointWriter>();
    auto tseries   = std::make_shared<TimeSeriesManager>();
    auto scheduler = std::make_shared<TaskScheduler>();
    double timestepSize = 0.0;

    auto bd = std::shared_ptr<BlockDecomposition>();
    auto bc = std::shared_ptr<BoundaryCondition> (new ThermalConvectionBoundaryCondition);
    auto rs = std::shared_ptr<RiemannSolver> (new HllcNewtonianHydroRiemannSolver);
    auto fs = std::shared_ptr<IntercellFluxScheme> (new MethodOfLinesPlm);
    auto cl = std::make_shared<NewtonianHydro>();


    // Set up grid shape. In 1D it's z-only. In 2D it's the x-z plane.
    // ------------------------------------------------------------------------
    auto mg = std::shared_ptr<MeshGeometry> (new SphericalMeshGeometry);
    //auto mg = std::shared_ptr<MeshGeometry> (new CartesianMeshGeometry);
    auto bs = Shape {{ 2, int (user["Nt"]) == 1 ? 0 : 2, 0, 0, 0 }};
    mg->setCellsShape ({{ user["Nr"], user["Nt"], 1 }});
    mg->setLowerUpper ({{ 1.0, M_PI * 0.5 - M_PI / 12.0, 0}}, {{10.0, M_PI * 0.5 + M_PI / 12.0, 0.1}});

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
    ss->setDisableFieldCT (true);
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

        for (unsigned int n = 0; n < fieldNames.size(); ++n)
        {
            entry[fieldNames[n]] = volumeAveraged[n];
        }
        tseries->append (status, entry);
    };

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
