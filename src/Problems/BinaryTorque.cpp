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




// Indexes to primitive quanitites P
#define RHO 0
#define V11 1
#define V22 2
#define V33 3
#define PRE 4
#define RAD 5

// Indexes to conserved quanitites U
#define DDD 0
#define S11 1
#define S22 2
#define S33 3
#define NRG 4




// ============================================================================
static const double SofteningRadius2 =  0.1;  // r0^2, where Fg = G M m / (r^2 + r0^2)
static const double BetaBuffer       = 0.01;  // Orbital periods over which to relax to IC in outer buffer
static const double BufferRadius     = 10.0;  // Radius beyond which solution is driven toward IC
static const double ViscousAlpha     = 5e-2;  // Alpha viscosity parameter
static const double MachNumber       = 20.0;  // 1/M = h/r
static const double GM               = 1.0;   // Newton G * mass
static const int    NumHoles         = 2;     // Number of Black holes 1 or 2
static const double aBin             = 1.0;
static const double OmegaBin         = 1.0;
static SimulationStatus status;




// ============================================================================
static double GravitationalPotential (double x, double y, double t)
{
    const double phi1 = 0.0;
    const double phi2 = M_PI;

    if (NumHoles == 1)
    {
        const double r = std::sqrt (x * x + y * y);
        return -GM / std::sqrt (r*r + SofteningRadius2);
    } 
    else if (NumHoles == 2)
    {
        const double x1 = 0.5 * aBin * std::cos(OmegaBin * t + phi1);
        const double y1 = 0.5 * aBin * std::sin(OmegaBin * t + phi1);
        const double x2 = 0.5 * aBin * std::cos(OmegaBin * t + phi2);
        const double y2 = 0.5 * aBin * std::sin(OmegaBin * t + phi2);
        const double r1s = std::sqrt((x-x1)*(x-x1) + (y-y1)*(y-y1) + SofteningRadius2);
        const double r2s = std::sqrt((x-x2)*(x-x2) + (y-y2)*(y-y2) + SofteningRadius2);
        return -(GM / r1s + GM / r2s);
    }
    else
    {
       assert (false);
    }
}

static std::array<double, 2> GravitationalAcceleration (double x, double y, double t)
{
    const double phi1 = 0.0;
    const double phi2 = M_PI;

    if (NumHoles == 1)
    {
        const double r = std::sqrt (x * x + y * y);
        const double xhat  = x / r;
        const double yhat  = y / r;
        return {-xhat * GM / (r * r + SofteningRadius2), -yhat * GM / (r * r + SofteningRadius2)};
    } 
    else if (NumHoles == 2)
    {
        const double x1 = 0.5 * aBin * std::cos(OmegaBin * t + phi1);
        const double y1 = 0.5 * aBin * std::sin(OmegaBin * t + phi1);
        const double x2 = 0.5 * aBin * std::cos(OmegaBin * t + phi2);
        const double y2 = 0.5 * aBin * std::sin(OmegaBin * t + phi2);
        const double r1 = std::sqrt((x-x1) * (x-x1) + (y-y1) * (y-y1));
        const double r2 = std::sqrt((x-x2) * (x-x2) + (y-y2) * (y-y2));
        const double r1s = std::sqrt((x-x1) * (x-x1) + (y-y1) * (y-y1) + SofteningRadius2);
        const double r2s = std::sqrt((x-x2) * (x-x2) + (y-y2) * (y-y2) + SofteningRadius2);
        const double xhat1 = x1 / r1;
        const double yhat1 = y1 / r1;
        const double xhat2 = x2 / r2;
        const double yhat2 = y2 / r2;
        const std::array<double, 2> F1 = {-xhat1 * GM / (r1s * r1s), -yhat1 * GM / (r1s * r1s)};
        const std::array<double, 2> F2 = {-xhat2 * GM / (r2s * r2s), -yhat2 * GM / (r2s * r2s)};
        return {F1[0] + F2[0], F1[1] + F2[1]};
    }
    else
    {
       assert (false);
    }
}




// ============================================================================
class ThinDiskNewtonianHydro : public ConservationLaw
{
public:
    ThinDiskNewtonianHydro();
    State fromConserved (const Request& request, const double* U) const override;
    State fromPrimitive (const Request& request, const double* P) const override;
    int getNumConserved() const override;
    int getIndexFor (VariableType type) const override;
    std::string getPrimitiveName (int fieldIndex) const override;
    std::vector<double> makeDiagnostics (const State& state) const override;
    std::vector<std::string> getDiagnosticNames() const override;
    double getSoundSpeedSquared (const double *P) const;
};




// ============================================================================
ThinDiskNewtonianHydro::ThinDiskNewtonianHydro()
{
}

ConservationLaw::State ThinDiskNewtonianHydro::fromConserved (const Request& request, const double* U) const
{
    double P[6];

    P[RHO] = U[DDD];
    P[PRE] = U[DDD] * getSoundSpeedSquared (P);
    P[V11] = U[S11] / U[DDD];
    P[V22] = U[S22] / U[DDD];
    P[V33] = U[S33] / U[DDD];
    P[RAD] = U[RAD];

    return fromPrimitive (request, P);
}

ConservationLaw::State ThinDiskNewtonianHydro::fromPrimitive (const Request& request, const double* P) const
{
    const auto dAA = request.areaElement;
    const double cs = std::sqrt (getSoundSpeedSquared (P));
    const double vn = P[V11] * dAA[0] + P[V22] * dAA[1] + P[V33] * dAA[2];

    auto S = State();
    S.numFields = 6;

    S.P[RHO] = P[RHO];
    S.P[V11] = P[V11];
    S.P[V22] = P[V22];
    S.P[V33] = P[V33];
    S.P[PRE] = P[PRE];
    S.P[RAD] = P[RAD];

    S.U[DDD] = P[RHO];
    S.U[S11] = P[RHO] * P[V11];
    S.U[S22] = P[RHO] * P[V22];
    S.U[S33] = P[RHO] * P[V33];
    S.U[NRG] = 0.0;
    S.U[RAD] = P[RAD];

    S.F[DDD] = vn * S.U[DDD];
    S.F[S11] = vn * S.U[S11] + P[PRE] * dAA[0];
    S.F[S22] = vn * S.U[S22] + P[PRE] * dAA[1];
    S.F[S33] = vn * S.U[S33] + P[PRE] * dAA[2];
    S.F[NRG] = 0.0;
    S.F[RAD] = 0.0;

    S.A[0] = vn - cs;
    S.A[1] = vn;
    S.A[2] = vn;
    S.A[3] = vn;
    S.A[4] = vn + cs;
    S.A[5] = 0.0;

    return S;
}

int ThinDiskNewtonianHydro::getNumConserved() const
{
    return 6;
}

int ThinDiskNewtonianHydro::getIndexFor (VariableType type) const
{
    switch (type)
    {
        case VariableType::density: return 0;
        case VariableType::velocity: return 1;
        case VariableType::pressure: return 4;
        default: return -1;
    }
}

std::string ThinDiskNewtonianHydro::getPrimitiveName (int fieldIndex) const
{
    switch (fieldIndex)
    {
        case RHO: return "density";
        case V11: return "velocity1";
        case V22: return "velocity2";
        case V33: return "velocity3";
        case PRE: return "pressure";
        case RAD: return "radius";
        default: return "";
    }
}

std::vector<double> ThinDiskNewtonianHydro::makeDiagnostics (const State& state) const
{
    auto D = std::vector<double> (5);
    D[0] = state.U[RHO];
    D[1] = state.U[S11];
    D[2] = state.U[S22];
    D[3] = state.U[S33];
    D[4] = state.U[NRG];
    return D;
}

std::vector<std::string> ThinDiskNewtonianHydro::getDiagnosticNames() const
{
    auto N = std::vector<std::string>(5);
    N[0] = "mass";
    N[1] = "momentum1";
    N[2] = "momentum2";
    N[3] = "momentum3";
    N[4] = "total_energy";
    return N;
}

double ThinDiskNewtonianHydro::getSoundSpeedSquared (const double* P) const
{
    return P[PRE] / P[DDD];
}




// ============================================================================
class ThinDiskBoundaryCondition : public BoundaryCondition
{
public:
    ThinDiskBoundaryCondition()
    {
    }

    void setConservationLaw (std::shared_ptr<ConservationLaw> cl) override
    {
    }

    void setMeshGeometry (std::shared_ptr<MeshGeometry> mg) override
    {
        meshGeometry = mg;
    }

    void apply (Cow::Array& A, MeshLocation location, MeshBoundary boundary, int axis, int numGuard) const override
    {
        resetPressureToIsothermal(A);
        outflow.apply (A, location, boundary, axis, numGuard);
    }

    bool isAxisPeriodic (int axis) override
    {
        return outflow.isAxisPeriodic (axis);
    }

    void resetPressureToIsothermal (Cow::Array& P) const
    {
        const double t = status.simulationTime;

        for (int i = 0; i < P.shape()[0]; ++i)
        {
            for (int j = 0; j < P.shape()[1]; ++j)
            {
                // WARNING: if other than 2 guard zones, do something!
                const auto X = meshGeometry->coordinateAtIndex (i - 2, j - 2, 0);
                const double r = std::sqrt (X[0] * X[0] + X[1] * X[1]);
                const double phi = GravitationalPotential (X[0], X[1], t);
                const double cs2 = std::pow (MachNumber, -2) * phi;
                const double sigma = P(i, j, 0, RHO);
                P(i, j, 0, PRE) = cs2 * sigma;
                P(i, j, 0, RAD) = r;
            }
        }
    }

    std::shared_ptr<MeshGeometry> meshGeometry;
    OutflowBoundaryCondition outflow;
};




// ============================================================================
static void computeViscousFluxes1D (const Array& P, double dx, Array& Fhat)
{
    for (int i = 0; i < P.shape()[0] - 1; ++i)
    {
        const double dgL = P(i + 0, 0, 0, RHO);
        const double dgR = P(i + 1, 0, 0, RHO);
        const double pgL = P(i + 0, 0, 0, PRE);
        const double pgR = P(i + 1, 0, 0, PRE);
        const double rcL = P(i + 0, 0, 0, RAD);
        const double rcR = P(i + 1, 0, 0, RAD);
        const double uyL = P(i + 0, 0, 0, V22);
        const double uyR = P(i + 1, 0, 0, V22);

        const double csL = std::sqrt (pgL / dgL);
        const double csR = std::sqrt (pgR / dgR);
        const double rc = 0.5 * (rcL + rcR);
        const double cs = 0.5 * (csL + csR);
        const double dg = 0.5 * (dgL + dgR);
        const double h0 = rc / MachNumber;
        const double nu = ViscousAlpha * cs * h0;
        const double mu = dg * nu;

        const double dxuy = (uyR - uyL) / dx;
        const double tauxy = mu * dxuy;

        Fhat(i + 1, 0, 0, 0, 0) -= 0.0;   // mass flux
        Fhat(i + 1, 0, 0, 1, 0) -= 0.0;   // px flux
        Fhat(i + 1, 0, 0, 2, 0) -= tauxy; // py flux
        Fhat(i + 1, 0, 0, 3, 0) -= 0.0;   // pz flux
        Fhat(i + 1, 0, 0, 4, 0) -= 0.0;   // E flux (neglected)
    }
}


static void computeViscousFluxes2D (const Array& P, double dx, double dy, Array& Fhat)
{
    for (int i = 0; i < P.shape()[0] - 1; ++i)
    {
        for (int j = 1; j < P.shape()[1] - 1; ++j)
        {
            const double dgL = P(i + 0, j, 0, RHO);
            const double dgR = P(i + 1, j, 0, RHO);
            const double pgL = P(i + 0, j, 0, PRE);
            const double pgR = P(i + 1, j, 0, PRE);
            const double rcL = P(i + 0, j, 0, RAD);
            const double rcR = P(i + 1, j, 0, RAD);
            const double uxL0 = P(i + 0, j - 1, 0, V11);
            const double uxL2 = P(i + 0, j + 0, 0, V11);
            const double uxR0 = P(i + 1, j + 0, 0, V11);
            const double uxR2 = P(i + 1, j + 1, 0, V11);
            const double uyL1 = P(i + 0, j - 1, 0, V22);
            const double uyR1 = P(i + 1, j + 1, 0, V22);

            const double csL = std::sqrt (pgL / dgL);
            const double csR = std::sqrt (pgR / dgR);
            const double rc = 0.5 * (rcL + rcR);
            const double cs = 0.5 * (csL + csR);
            const double dg = 0.5 * (dgL + dgR);
            const double h0 = rc / MachNumber;
            const double nu = ViscousAlpha * cs * h0;
            const double mu = dg * nu;

            const double dxuy = (uyR1 - uyL1) / dx;
            const double dyux = (uxR2 - uxR0 + uxL2 - uxL0) / (4 * dy);

            // To be checked...
            const double tauxx = 0.0;
            const double tauxy = mu * (dxuy + dyux);

            Fhat(i + 1, j, 0, 0, 0) -= 0.0;   // mass flux
            Fhat(i + 1, j, 0, 1, 0) -= tauxx; // px flux
            Fhat(i + 1, j, 0, 2, 0) -= tauxy; // py flux
            Fhat(i + 1, j, 0, 3, 0) -= 0.0;   // pz flux
            Fhat(i + 1, j, 0, 4, 0) -= 0.0;   // E flux (neglected)
        }
    }

    for (int i = 1; i < P.shape()[0] - 1; ++i)
    {
        for (int j = 0; j < P.shape()[1] - 1; ++j)
        {
            const double dgL = P(i, j + 0, 0, RHO);
            const double dgR = P(i, j + 1, 0, RHO);
            const double pgL = P(i, j + 0, 0, PRE);
            const double pgR = P(i, j + 1, 0, PRE);
            const double rcL = P(i, j + 0, 0, RAD);
            const double rcR = P(i, j + 1, 0, RAD);
            const double uyL0 = P(i - 1, j + 0, 0, V22);
            const double uyL2 = P(i + 0, j + 0, 0, V22);
            const double uyR0 = P(i + 0, j + 1, 0, V22);
            const double uyR2 = P(i + 1, j + 1, 0, V22);
            const double uxL1 = P(i - 1, j + 0, 0, V11);
            const double uxR1 = P(i + 1, j + 1, 0, V11);

            const double csL = std::sqrt (pgL / dgL);
            const double csR = std::sqrt (pgR / dgR);
            const double rc = 0.5 * (rcL + rcR);
            const double cs = 0.5 * (csL + csR);
            const double dg = 0.5 * (dgL + dgR);
            const double h0 = rc / MachNumber;
            const double nu = ViscousAlpha * cs * h0;
            const double mu = dg * nu;

            const double dyux = (uxR1 - uxL1) / dy;
            const double dxuy = (uyR2 - uyR0 + uyL2 - uyL0) / (4 * dx);

            // To be checked...
            const double tauyx = mu * (dyux + dxuy);
            const double tauyy = 0.0;

            Fhat(i, j + 1, 0, 0, 1) -= 0.0;   // mass flux
            Fhat(i, j + 1, 0, 1, 1) -= tauyx; // px flux
            Fhat(i, j + 1, 0, 2, 1) -= tauyy; // py flux
            Fhat(i, j + 1, 0, 3, 1) -= 0.0;   // pz flux
            Fhat(i, j + 1, 0, 4, 1) -= 0.0;   // E flux (neglected)
        }
    }
}

static auto makeViscousFluxCorrection (double dx, double dy)
{
    return [dx,dy] (const Cow::Array& P, Cow::Array& F)
    {
        if (P.shape()[1] == 1)
        {
            computeViscousFluxes1D (P, dx, F);
        }
        else
        {
            computeViscousFluxes2D (P, dx, dy, F);
        }
    };
}




// ============================================================================
int BinaryTorque::run (int argc, const char* argv[])
{
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
    user["theta"]   = 0.0; // angle of shear layer (in degrees)

    auto cl = std::make_shared<ThinDiskNewtonianHydro>();


    // Initial data function for testing
    // ------------------------------------------------------------------------
    auto initialDataLinearShear = [&] (double x, double y, double z) -> std::vector<double>
    {
        const double theta = double (user["theta"]) * M_PI / 180.0;
        const double nx = std::cos (theta);
        const double ny = std::sin (theta);
        const double r = x * nx + y * ny;
        const double nu = 0.1;
        const double t0 = 1.0 / nu;
        const double ur = std::exp (-r * r / (4 * nu * t0));
        return std::vector<double> {1.0, -ur * ny, ur * nx, 0.0, 1.0, 0.0};
    };


    auto initialData = [&] (double x, double y, double z) -> std::vector<double>
    {
        const double t = status.simulationTime;
        const double r2 = x * x + y * y;
        const double r  = std::sqrt (r2);
        const double rho = 0.1 + std::exp (-std::pow (r - 4, 2));
        const double vq = std::sqrt (-GravitationalPotential (x, y, t));
        const double pre = std::pow (vq / MachNumber, 2.0) * rho; // cs^2 * rho
        const double qhX = -y / r;     
        const double qhY =  x / r;
        return std::vector<double> {rho, vq * qhX, vq * qhY, 0, pre, 0.0};
    };


    // Source terms
    // ------------------------------------------------------------------------
    auto sourceTermsFunction = [&] (double x, double y, double z, StateArray P)
    {
        const double t = status.simulationTime;
        const double dg = P[0];
        const double vx = P[1];
        const double vy = P[2];
        const double r = std::sqrt (x * x + y * y);
        const auto ag = GravitationalAcceleration (x, y, t);

        auto S = StateArray();

        // Gravitational source term
        // ====================================================================
        S[DDD] = 0.0;
        S[S11] = dg * ag[0];
        S[S22] = dg * ag[1];
        S[S33] = 0.0;
        S[NRG] = dg * (ag[0] * vx + ag[1] * vy);
        S[RAD] = 0.0;

        ConservationLaw::Request request;
        request.getPrimitive = true;
        request.getConserved = true;
        request.getFluxes = false;
        request.getEigenvalues = false;
        request.areaElement[0] = 1.0;
        request.areaElement[1] = 0.0;
        request.areaElement[2] = 0.0;

        auto state0 = cl->fromPrimitive (request, &initialData (x, y, z)[0]);
        auto state1 = cl->fromPrimitive (request, &P[0]);

        const double vk   = std::sqrt (-GravitationalPotential (x, y, t));
        const double torb = 2.0 * M_PI * r / vk;
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
    // auto rs = std::shared_ptr<RiemannSolver> (new HlleRiemannSolver);
    auto fs = std::shared_ptr<IntercellFluxScheme> (new MethodOfLinesPlm);
    auto bc = std::shared_ptr<BoundaryCondition>(new ThinDiskBoundaryCondition);

    // Set up grid shape.
    // ------------------------------------------------------------------------
    int N = int (user["N"]);
    auto cs = Shape {{ N, N, 1, 1, 1 }};
    auto bs = Shape {{ 2, 2, 0, 0, 0 }};
    // auto cs = Shape {{ N, 1, 1, 1, 1 }};
    // auto bs = Shape {{ 2, 0, 0, 0, 0 }};

    mg->setCellsShape (cs);
    mg->setLowerUpper ({{-10, -10, 0.0}}, {{10, 10, 1.0}});
    // mg->setLowerUpper ({{-10, 0.0, 0.0}}, {{10, 1.0, 1.0}});

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


    double dx = mg->cellLength (0, 0, 0, 0);
    double dy = mg->cellLength (0, 0, 0, 1);

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
    ss->setViscousFluxFunction (makeViscousFluxCorrection (dx, dy));

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
        // md->assignPrimitive (mo->generate (initialDataLinearShear, MeshLocation::cell));
    }
    md->applyBoundaryCondition (*bc);
    taskRecomputeDt (status, 0);

    return maraMainLoop (status, timestep, condition, advance, *scheduler, *logger);  
}
