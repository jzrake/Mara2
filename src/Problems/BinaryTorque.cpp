#include <cmath>
#include <iomanip>
#include "BinaryTorque.hpp"
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




// ============================================================================
static const double SofteningEpsilon =  0.1;  // r0^2, where Fg = G M m / (r^2 + r0^2)
static const double BetaBuffer       = 0.01;  // Orbital periods over which to relax to IC in outer buffer
static const double BufferRadius     = 10.0;  // Radius beyond which solution is driven toward IC
static const double ViscousAlpha     =  0.1;  // Alpha viscosity parameter




// ============================================================================
int BinaryTorque::run (int argc, const char* argv[])
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
    user["plm"]     = 1.8;
    user["N"]       = 128;

    auto cl = std::make_shared<NewtonianHydro>();


    auto initialData = [&] (double x, double y, double z) -> std::vector<double>
    {
        const double GM = 1.0;
        const double r2 = x * x + y * y;
        const double r  = std::sqrt (r2);
        const double rho = 1;
        const double pre = 1;
        const double vq = std::sqrt (GM * r / (r2 + SofteningEpsilon));
        const double qhX = -y / r;
        const double qhY =  x / r;

        return std::vector<double> {rho, vq * qhX, vq * qhY, 0, pre};
    };


    // Source terms
    // ------------------------------------------------------------------------
    auto sourceTermsFunction = [&] (double x, double y, double z, StateArray P)
    {
        const double GM = 1.0;
        const double dg = P[0];
        const double vx = P[1];
        const double vy = P[2];
        const double r2 = x * x + y * y;
        const double r  = std::sqrt (r2);
        const double xh = x / r;
        const double yh = y / r;
        const double vr = vx * xh + vy * yh;
        const double fg = -GM * dg / (r2 + SofteningEpsilon);

        auto S = StateArray();

        // Gravitational source term
        // ====================================================================
        S[0] = 0.0;
        S[1] = fg * xh;
        S[2] = fg * yh;
        S[3] = 0.0;
        S[4] = fg * vr;

        NewtonianHydro::Request request;
        request.getPrimitive = true;
        request.getConserved = true;
        request.getFluxes = false;
        request.getEigenvalues = false;
        request.areaElement[0] = 1.0;
        request.areaElement[1] = 0.0;
        request.areaElement[2] = 0.0;

        auto state0 = cl->fromPrimitive (request, &initialData (x, y, z)[0]);
        auto state1 = cl->fromPrimitive (request, &P[0]);

        const double torb = 2.0 * M_PI / std::sqrt (GM / r / (r2 + SofteningEpsilon));
        const double env  = 2.0 / (1 + std::tanh (r - BufferRadius));
        const double tau  = BetaBuffer * env * torb;

        // Outer buffer
        // ====================================================================
        if (r > 0.5 * BufferRadius)
        {
            S[0] -= (state1.U[0] - state0.U[0]) / tau;
            S[1] -= (state1.U[1] - state0.U[1]) / tau;
            S[2] -= (state1.U[2] - state0.U[2]) / tau;
            S[3] -= (state1.U[3] - state0.U[3]) / tau;
            S[4] -= (state1.U[4] - state0.U[4]) / tau;
        }

        // Sink radius
        // ====================================================================
        if (r < 1.0)
        {
            S[0] -= state1.U[0] / (torb / ViscousAlpha);
            S[1] -= state1.U[1] / (torb / ViscousAlpha);
            S[2] -= state1.U[2] / (torb / ViscousAlpha);
            S[3] -= state1.U[3] / (torb / ViscousAlpha);
            S[4] -= state1.U[4] / (torb / ViscousAlpha);
        }

        return S;
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
    auto mg = std::shared_ptr<MeshGeometry> (new CartesianMeshGeometry);
    auto rs = std::shared_ptr<RiemannSolver> (new HllcNewtonianHydroRiemannSolver);
    auto fs = std::shared_ptr<IntercellFluxScheme> (new MethodOfLinesPlm);
    auto bc = std::shared_ptr<BoundaryCondition>(new PeriodicBoundaryCondition);

    // Set up grid shape.
    // ------------------------------------------------------------------------
    int N = int (user["N"]);
    auto cs = Shape {{ N, N, 1, 1, 1 }};
    auto bs = Shape {{ 2, 2, 0, 0, 0 }};

    mg->setCellsShape (cs);
    mg->setLowerUpper ({{-10, -10, 0.0}}, {{10, 10, 1.0}});

    if (! user["serial"])
    {
        logger->setLogToNullUnless (MpiCommunicator::world().isThisMaster());
        bd = std::make_shared<BlockDecomposition> (mg, *logger);
        mg = bd->decompose();
        bc = bd->createBoundaryCondition (bc);
    }
    bc->setMeshGeometry (mg);
    bc->setConservationLaw (cl);

    tseries->setLogger (logger);
    scheduler->setLogger (logger);

    auto ss = std::make_shared<MethodOfLinesTVD>();
    auto mo = std::make_shared<MeshOperator>();
    auto fo = std::make_shared<FieldOperator>();
    auto md = std::make_shared<MeshData> (mg->cellsShape(), bs, cl->getNumConserved());

    cl->setGammaLawIndex (4. / 3);
    // cl->setPressureFloor (1e-2);
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

    scheduler->schedule (taskCheckpoint, TaskScheduler::Recurrence (user["cpi"]).logarithmic(), "checkpoint");
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

    return maraMainLoop (status, timestep, condition, advance, *scheduler, *logger);  
}
