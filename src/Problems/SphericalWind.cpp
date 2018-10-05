#include <cmath>
#include "SphericalWind.hpp"
#include "../Mara.hpp"
#include "../BoundaryConditions.hpp"
#include "../BlockDecomposition.hpp"
#include "../SphericalMeshGeometry.hpp"
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
class JetLaunchingBoundaryCondition : public BoundaryCondition
{
public:
    static double structureExponent;
    static double jetOpeningAngle;
    static double jetGamma;
    static double jetEnthalpy;

    void setConservationLaw (std::shared_ptr<ConservationLaw> cl) override
    {
        reflecting.setConservationLaw (cl);
    }

    void setMeshGeometry (std::shared_ptr<MeshGeometry> mg) override
    {
        meshGeometry = mg;
    }

    void apply (Cow::Array& A, MeshLocation location, MeshBoundary boundary, int axis, int numGuard) const override
    {
        if (axis == 1)
        {
            return reflecting.apply (A, location, boundary, axis, numGuard);
        }
        if (boundary == MeshBoundary::left)
        {
            return applyInflowAtInnerBoundary (A, numGuard);
        }
        return outflow.apply (A, location, boundary, axis, numGuard);
    }

    void applyInflowAtInnerBoundary (Cow::Array& P, int numGuard) const
    {
        for (int i = 0; i < numGuard; ++i)
        {
            for (int j = numGuard; j < P.size(1) - numGuard; ++j)
            {
                const double q = meshGeometry->coordinateAtIndex (i - numGuard, j - numGuard, 0)[1];

                const double alpha  = structureExponent;
                const double qjet   = jetOpeningAngle;
                const double w0     = jetEnthalpy;
                const double fjet   = std::exp (-std::pow (q / qjet, alpha));
                const double gamma0 = 1.0 + jetGamma * fjet;

                const double vr = std::sqrt (1.0 - 1.0 / (gamma0 * gamma0));
                const double rho = 1.0;
                const double pre = rho * w0 / 4;

                P(i, j, 0, 0) = rho;
                P(i, j, 0, 1) = vr;
                P(i, j, 0, 2) = 0.0;
                P(i, j, 0, 3) = 0.0;
                P(i, j, 0, 4) = pre;
                P(i, j, 0, 5) = 0.0;
                P(i, j, 0, 6) = 0.0;
                P(i, j, 0, 7) = 0.0;
            }
        }
    }

    bool isAxisPeriodic (int axis) override
    {
        switch (axis)
        {
            case 0: return false;
            case 1: return false;
            case 2: return true;
            default: throw std::logic_error ("SphericalWind");
        }
    }
    std::shared_ptr<MeshGeometry> meshGeometry;
    OutflowBoundaryCondition outflow;
    ReflectingBoundaryCondition reflecting;
};

double JetLaunchingBoundaryCondition::structureExponent = 1.5;
double JetLaunchingBoundaryCondition::jetOpeningAngle   = 0.2;
double JetLaunchingBoundaryCondition::jetGamma          = 5.0;
double JetLaunchingBoundaryCondition::jetEnthalpy       = 5.0;




// ============================================================================
int SphericalWind::run (int argc, const char* argv[])
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
    user["plm"]     = 1.5;
    user["Nr"]      = 128;
    user["Nt"]      = 1;
    user["router"]  = 1e3;
    user["qtor"]    = 0.0;
    user["dtor"]    = 1e2;
    user["d0"]      = 1.0;
    user["dind"]    = 3.0;
    user["u0"]      = 1e-1;
    user["peng"]    = 20;
    user["alpha"]   = 1.5;
    user["theta0"]  = 0.2;
    user["gamma0"]  = 5.0;
    user["w0"]      = 5.0;



    // Geometrical source terms
    // ------------------------------------------------------------------------
    auto sourceTermsFunction = [&] (double r, double q, double p, StateArray P)
    {
        const double gm = 4. / 3;
        const double dg = P[0];
        const double vr = P[1];
        const double vq = P[2];
        const double vp = P[3];
        const double pg = P[4];
        const double eg = pg / dg / (gm - 1);
        const double v2 = vr * vr + vq * vq + vp * vp;
        const double W2 = 1.0 / (1.0 - v2);
        const double hg = 1.0 + eg + pg / dg;
        const double rhohW2 = dg * hg * W2;

        auto S = StateArray();

        S[0] = 0.0;
        S[1] = (2 * pg + rhohW2 * (vq * vq + vp * vp)) / r;
        S[2] = (pg * cot(q) + rhohW2 * (vp * vp * cot(q) - vr * vq)) / r;
        S[3] = -rhohW2 * vp * (vr + vq * cot(q)) / r;
        S[4] = 0.0;

        S[5] = 0.0;
        S[6] = 0.0;
        S[7] = 0.0;

        return S;
    };

    auto initialData = [&] (double r, double q, double p) -> std::vector<double>
    {
        double d0   = user["d0"];
        double dind = user["dind"];
        double rho  = d0 * std::pow (r, -dind);
        double pre  = std::pow (rho, 4. / 3) * 1e-2;
        double ur   = r * double (user["u0"]);
        double vr   = ur / std::sqrt (1 + ur * ur);
        return std::vector<double> {rho, vr, 0, 0, pre, 0, 0, 0};
    };

    auto clineUser = Variant::fromCommandLine (argc, argv);
    auto chkptUser = Variant::NamedValues();

    if (clineUser.find ("restart") != clineUser.end())
    {
        chkptUser = H5::File (clineUser["restart"], "r").getGroup ("user").readNamedValues();
    }
    Variant::update (user, chkptUser);
    Variant::update (user, clineUser);

    JetLaunchingBoundaryCondition::structureExponent = user["alpha"];
    JetLaunchingBoundaryCondition::jetOpeningAngle   = user["theta0"];
    JetLaunchingBoundaryCondition::jetGamma          = user["gamma0"];
    JetLaunchingBoundaryCondition::jetEnthalpy       = user["w0"];

    auto logger    = std::make_shared<Logger>();
    auto writer    = std::make_shared<CheckpointWriter>();
    auto tseries   = std::make_shared<TimeSeriesManager>();
    auto scheduler = std::make_shared<TaskScheduler>();

    double timestepSize = 0.0;

    auto bd = std::shared_ptr<BlockDecomposition>();
    auto mg = std::shared_ptr<MeshGeometry> (new SphericalMeshGeometry);
    auto rs = std::shared_ptr<RiemannSolver> (new HlleRiemannSolver);
    auto fs = std::shared_ptr<IntercellFluxScheme> (new MethodOfLinesPlm);
    auto bc = std::shared_ptr<BoundaryCondition>(new JetLaunchingBoundaryCondition);
    auto cl = std::make_shared<RelativisticMHD>();

    // Set up grid shape.
    // ------------------------------------------------------------------------
    int Nr = int (user["Nr"]);
    int Nt = int (user["Nt"]);
    auto cs = Shape {{ Nr, Nt, 1, 1, 1 }};
    auto bs = Shape {{ 2, Nt == 1 ? 0 : 2, 0, 0, 0 }};

    mg->setCellsShape (cs);

    if (Nt == 1)
        mg->setLowerUpper ({{1.0, M_PI_2 - 0.1, 0.0}}, {{user["router"], M_PI_2 + 0.1, 0.1}});
    else
        mg->setLowerUpper ({{1.0, 0.0, 0.0}}, {{user["router"], M_PI_2, 0.1}});

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
    cl->setPressureFloor (1e-2);
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

    logger->log() << std::endl << user << std::endl;
    return maraMainLoop (status, timestep, condition, advance, *scheduler, *logger);  
}






// Initial conditions function for torus-cloud:

// auto initialData = [&] (double r, double q, double p) -> std::vector<double>
// {
//     double qtor = user["qtor"];
//     double dtor = user["dtor"];
//     double peng = user["peng"];
//     double rhor = std::pow (r, -double (user["dind"]));
//     double rhoq = rhor * (qtor == 0.0 ? 1.0 : std::exp (-r) * (1 + dtor * std::exp (-std::pow ((q - M_PI_2) / qtor, 2))));
//     double pre = std::pow (rhor, 4. / 3) * 1e-3 + peng * std::exp (-r);
//     double ur = r * double (user["u0"]);
//     double vr = ur / std::sqrt (1 + ur * ur);
//     return std::vector<double> {rhoq, vr, 0, 0, pre, 0, 0, 0};
// };