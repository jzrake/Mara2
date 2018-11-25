#include <cmath>
#include <cassert>
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
#define HUGE_TIME_SCALE 1e20
#define DENSITY_FLOOR 1e-6

using namespace Cow;




// Indexes to primitive quanitites P
#define RHO 0
#define V11 1
#define V22 2
#define V33 3
#define PRE 4
#define XXX 5
#define YYY 6

// Indexes to conserved quanitites U
#define DDD 0
#define S11 1
#define S22 2
#define S33 3
#define NRG 4




// ============================================================================
static double SofteningRadius  =  0.2;  // r0^2, where Fg = G M m / (r^2 + r0^2)
static double BetaBuffer       =  1.0;  // Orbital periods over which to relax to IC in outer buffer
static double BufferRadius     = 10.0;  // Radius beyond which solution is driven toward IC
static double ViscousAlpha     = 1e-1;  // Alpha viscosity parameter
static double MachNumber       = 10.0;  // 1/M = h/r
static int    NumHoles         = 1;     // Number of Black holes 1 or 2
static double GM               = 1.0;   // Newton G * mass
static double aBin             = 1.0;   // Binary separation
static double SinkRadius       = 0.1;   // Sink radius

// Rotating frame of reference

static double omegaFrame = 0.0;

static std::string InitialDataString = "Ring";

static SimulationStatus status;




// ============================================================================
struct SinkGeometry
{
    double x;
    double y;
    double r;
    double xhat;
    double yhat;
};

static std::vector<SinkGeometry> SinkGeometries (double x, double y, double t)
{
    if (NumHoles == 1)
    {
        SinkGeometry g;
        g.x = x;
        g.y = y;
        g.r = std::sqrt (x * x + y * y);
        g.xhat = g.x / g.r;
        g.yhat = g.y / g.r;
        return {g};
    }
    if (NumHoles == 2)
    {
        const double omegaBin = aBin == 0.0 ? 0.0 : std::sqrt (2.0 * GM / (aBin * aBin * aBin));
        SinkGeometry g1, g2;
        g1.x = 0.5 * aBin * std::cos ((omegaBin - omegaFrame) * t);
        g1.y = 0.5 * aBin * std::sin ((omegaBin - omegaFrame) * t);
        g2.x = 0.5 * aBin * std::cos ((omegaBin - omegaFrame) * t + M_PI);
        g2.y = 0.5 * aBin * std::sin ((omegaBin - omegaFrame) * t + M_PI);
        g1.r = std::sqrt ((x - g1.x) * (x - g1.x) + (y - g1.y) * (y - g1.y));
        g2.r = std::sqrt ((x - g2.x) * (x - g2.x) + (y - g2.y) * (y - g2.y));
        g1.xhat = (x - g1.x) / g1.r;
        g1.yhat = (y - g1.y) / g1.r;
        g2.xhat = (x - g2.x) / g2.r;
        g2.yhat = (y - g2.y) / g2.r;
        return {g1, g2};
    }
    return {};
}

static double GravitationalPotential (double x, double y, double t)
{
    const double rs = SofteningRadius;
    double phi = 0.0;

    for (const auto& g : SinkGeometries (x, y, t))
    {
        phi += -GM / (g.r + rs);
    }
    return phi;
}

static std::array<double, 2> GravitationalAcceleration (double x, double y, double t)
{
    const double rs = SofteningRadius;
    std::array<double, 2> a = {{0.0, 0.0}};

    for (const auto& g : SinkGeometries (x, y, t))
    {
        a[0] += -g.xhat * GM / std::pow (g.r + rs, 2);
        a[1] += -g.yhat * GM / std::pow (g.r + rs, 2);
    }
    return a;
}

static double SinkKernel (double r)
{
    const double viscousTimeScale = 2 * M_PI / std::sqrt (SinkRadius * SinkRadius * SinkRadius / GM) / ViscousAlpha;
    return r < SinkRadius ? viscousTimeScale : HUGE_TIME_SCALE;
}

static double SinkTime (double x, double y, double t)
{
    double tmin = HUGE_TIME_SCALE;

    for (const auto& g : SinkGeometries (x, y, t))
    {
        if (SinkKernel (g.r) < tmin)
        {
            tmin = SinkKernel (g.r);
        }
    }
    return tmin;
}

static double bufferZoneProfile (double r)
{
    return 1 + std::tanh (r - BufferRadius);
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
    double P[7] = {0, 0, 0, 0, 0, 0, 0};

    const double t = status.simulationTime;
    const double x = U[XXX];
    const double y = U[YYY];
    const double phi = GravitationalPotential (x, y, t);
    const double cs2 = std::pow (MachNumber, -2) * -phi;

    P[RHO] = std::max (DENSITY_FLOOR, U[DDD]);
    P[PRE] = P[RHO] * cs2;
    P[V11] = U[S11] / P[RHO];
    P[V22] = U[S22] / P[RHO];
    P[V33] = U[S33] / P[RHO];
    P[XXX] = U[XXX];
    P[YYY] = U[YYY];

    if (P[PRE] < 0.0 || P[RHO] < 0.0)
    {
        auto S = State();
        S.numFields = getNumConserved();

        for (int q = 0; q < getNumConserved(); ++q)
        {
            S.U[q] = U[q];
            S.P[q] = P[q];
        }
        throw ConservationLaw::StateFailure (S);
    }
    return fromPrimitive (request, P);
}

ConservationLaw::State ThinDiskNewtonianHydro::fromPrimitive (const Request& request, const double* P) const
{
    const auto dAA = request.areaElement;
    const double cs = std::sqrt (getSoundSpeedSquared (P));
    const double vn = P[V11] * dAA[0] + P[V22] * dAA[1] + P[V33] * dAA[2];

    auto S = State();
    S.numFields = 7;

    S.P[RHO] = P[RHO];
    S.P[V11] = P[V11];
    S.P[V22] = P[V22];
    S.P[V33] = P[V33];
    S.P[PRE] = P[PRE];
    S.P[XXX] = P[XXX];
    S.P[YYY] = P[YYY];

    S.U[DDD] = P[RHO];
    S.U[S11] = P[RHO] * P[V11];
    S.U[S22] = P[RHO] * P[V22];
    S.U[S33] = P[RHO] * P[V33];
    S.U[NRG] = P[PRE];
    S.U[XXX] = P[XXX];
    S.U[YYY] = P[YYY];

    S.F[DDD] = vn * S.U[DDD];
    S.F[S11] = vn * S.U[S11] + P[PRE] * dAA[0];
    S.F[S22] = vn * S.U[S22] + P[PRE] * dAA[1];
    S.F[S33] = vn * S.U[S33] + P[PRE] * dAA[2];
    S.F[NRG] = 0.0;
    S.F[XXX] = 0.0;
    S.F[YYY] = 0.0;

    S.A[0] = vn - cs;
    S.A[1] = vn;
    S.A[2] = vn;
    S.A[3] = vn;
    S.A[4] = vn + cs;
    S.A[XXX] = 0.0;
    S.A[YYY] = 0.0;

    return S;
}

int ThinDiskNewtonianHydro::getNumConserved() const
{
    return 7;
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
        case XXX: return "x";
        case YYY: return "y";
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
    return P[PRE] / P[RHO];
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
        return false;
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
                const double phi = GravitationalPotential (X[0], X[1], t);
                const double cs2 = std::pow (MachNumber, -2) * -phi;
                const double sigma = P(i, j, 0, RHO);
                P(i, j, 0, PRE) = cs2 * sigma;
                P(i, j, 0, XXX) = X[0];
                P(i, j, 0, YYY) = X[1];
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
        const double xcL = P(i + 0, 0, 0, XXX);
        const double xcR = P(i + 1, 0, 0, XXX);
        const double ycL = P(i + 0, 0, 0, YYY);
        const double ycR = P(i + 1, 0, 0, YYY);
        const double uyL = P(i + 0, 0, 0, V22);
        const double uyR = P(i + 1, 0, 0, V22);
        const double rcL = std::sqrt (xcL * xcL + ycL * ycL);
        const double rcR = std::sqrt (xcR * xcR + ycR * ycR);

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
            const double xcL = P(i + 0, j, 0, XXX);
            const double xcR = P(i + 1, j, 0, XXX);
            const double ycL = P(i + 0, j, 0, YYY);
            const double ycR = P(i + 1, j, 0, YYY);
            const double uxL0 = P(i + 0, j - 1, 0, V11);
            const double uxL2 = P(i + 0, j + 0, 0, V11);
            const double uxR0 = P(i + 1, j + 0, 0, V11);
            const double uxR2 = P(i + 1, j + 1, 0, V11);
            const double uyL1 = P(i + 0, j - 1, 0, V22);
            const double uyR1 = P(i + 1, j + 1, 0, V22);
            const double rcL = std::sqrt (xcL * xcL + ycL * ycL);
            const double rcR = std::sqrt (xcR * xcR + ycR * ycR);

            const double csL = std::sqrt (pgL / dgL);
            const double csR = std::sqrt (pgR / dgR);
            const double rc = 0.5 * (rcL + rcR);
            const double cs = 0.5 * (csL + csR);
            const double dg = 0.5 * (dgL + dgR);
            const double h0 = rc / MachNumber;
            const double nu = ViscousAlpha * cs * h0;// * (rc < 8.0 ? 1.0 : 0.0);
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
            const double xcL = P(i, j + 0, 0, XXX);
            const double xcR = P(i, j + 1, 0, XXX);
            const double ycL = P(i, j + 0, 0, YYY);
            const double ycR = P(i, j + 1, 0, YYY);
            const double uyL0 = P(i - 1, j + 0, 0, V22);
            const double uyL2 = P(i + 0, j + 0, 0, V22);
            const double uyR0 = P(i + 0, j + 1, 0, V22);
            const double uyR2 = P(i + 1, j + 1, 0, V22);
            const double uxL1 = P(i - 1, j + 0, 0, V11);
            const double uxR1 = P(i + 1, j + 1, 0, V11);
            const double rcL = std::sqrt (xcL * xcL + ycL * ycL);
            const double rcR = std::sqrt (xcR * xcR + ycR * ycR);

            const double csL = std::sqrt (pgL / dgL);
            const double csR = std::sqrt (pgR / dgR);
            const double rc = 0.5 * (rcL + rcR);
            const double cs = 0.5 * (csL + csR);
            const double dg = 0.5 * (dgL + dgR);
            const double h0 = rc / MachNumber;
            const double nu = ViscousAlpha * cs * h0 * (rc < 8.0 ? 1.0 : 0.0);
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
    user["cfl"]     = 0.3;
    user["plm"]     = 1.5;
    user["N"]       = 128;

    user["SofteningRadius"]  = SofteningRadius;
    user["BetaBuffer"]       = BetaBuffer;
    user["BufferRadius"]     = BufferRadius;
    user["ViscousAlpha"]     = ViscousAlpha;
    user["MachNumber"]       = MachNumber;
    user["NumHoles"]         = NumHoles;
    user["aBin"]             = aBin;
    user["SinkRadius"]       = SinkRadius;
    user["InitialData"]      = InitialDataString;

    auto cl = std::make_shared<ThinDiskNewtonianHydro>();
    auto initialData = std::function<std::vector<double>(double x, double y, double z)>();
    double timestepSize = 0.0;


    // Initial data function for testing
    // ------------------------------------------------------------------------
    auto initialDataLinearShear = [&] (double x, double y, double z) -> std::vector<double>
    {
        const double theta = 0.0;
        const double nx = std::cos (theta);
        const double ny = std::sin (theta);
        const double r = x * nx + y * ny;
        const double nu = 0.1;
        const double t0 = 1.0 / nu;
        const double ur = std::exp (-r * r / (4 * nu * t0));
        return std::vector<double> {1.0, -ur * ny, ur * nx, 0.0, 1.0, x, y};
    };

    auto initialDataRing = [&] (double x, double y, double z) -> std::vector<double>
    {
        const double t = status.simulationTime;
        const double r2 = x * x + y * y;
        const double r  = std::sqrt (r2);
        const double rho = 0.1 + std::exp (-std::pow (r - 4, 2));
        const double vq = std::sqrt (-GravitationalPotential (x, y, t));
        const double pre = std::pow (vq / MachNumber, 2.0) * rho; // cs^2 * rho
        const double qhX = -y / r;
        const double qhY =  x / r;
        return std::vector<double> {rho, vq * qhX, vq * qhY, 0, pre, x, y};
    };

    auto initialDataFarris14 = [&] (double x, double y, double z) -> std::vector<double>
    {
        // Initial conditions from Farris+ (2014) ApJ 783, 134
        // MISSING: radial drift velocity, pressure gradient for omega
        const double rs = 10.0 * aBin;
        const double delta = 3.0;
        const double xsi = 2.0;
        const double sigma0 = 1.0;
        const double r2 = x * x + y * y;
        const double r  = std::sqrt (r2);
        const double rho = sigma0 * std::pow ((r/rs), -delta) * std::exp (-std::pow (r/rs, -xsi));
        const double omegak2 = GM/(r*r2); // squared Keplerian for total M at center, is M total mass of binary?
        const double omega2  = omegak2 * std::pow (1.0 + (3.0/16.0)*aBin*aBin/r2, 2.0); //missing pressure gradient term
        const double vq = r * std::sqrt (omega2);
        const double pre = std::pow (vq / MachNumber, 2.0) * rho; // cs^2 * rho
        const double qhX = -y / r;
        const double qhY =  x / r;
        return std::vector<double> {rho, vq * qhX, vq * qhY, 0, pre, x, y};
    };

    auto initialDataTang17 = [&] (double x, double y, double z) -> std::vector<double>
    {
        // Initial conditions from Tang+ (2017) MNRAS 469, 4258
        // MISSING?: radial drift velocity, pressure gradient for omegak2
        // Check Yike/Chris Disco setup
        const double rs = SofteningRadius;
        const double r0 = 2.5 * aBin;
        const double xsi = 10.0;
        const double sigma0 = 1.0;
        const double r2 = x * x + y * y;
        const double r  = std::sqrt (r2);
        const double sigma = sigma0 * std::pow (r + rs, -0.5) * std::max (std::exp (-std::pow (r/r0, -xsi)), 1e-6);
        const double omegak2 = GM / std::pow (r + rs,  3.0) * NumHoles; // squared Keplerian for total M at center, is M total mass of binary?
        const double omega2  = omegak2;// * std::pow (1.0 + (3.0/16.0)*aBin*aBin/r2, 2.0); // missing pressure gradient term
        const double vq = r * std::sqrt (omega2);
        const double pre = std::pow (vq / MachNumber, 2.0) * sigma; // cs^2 * rho
        const double h0 = r / MachNumber;
        const double nu = ViscousAlpha * (vq / MachNumber) * h0; //ViscousAlpha * cs * h0
        const double vr = -(3.0/2.0) * nu / r; //radial drift velocity

        const double vx = vq * (-y / r) + vr * (x / r);
        const double vy = vq * ( x / r) + vr * (y / r);

        return std::vector<double> {sigma, vx, vy, 0, pre, x, y};
    };

   auto initialDataTang17Rot = [&] (double x, double y, double z) -> std::vector<double>
    {
        // Initial conditions from Tang+ (2017) MNRAS 469, 4258
        // MISSING?: radial drift velocity, pressure gradient for omegak2
        // Check Yike/Chris Disco setup
        const double rs = SofteningRadius;
        const double r0 = 2.5 * aBin;
        const double xsi = 10.0;
        const double sigma0 = 1.0;
        const double r2 = x * x + y * y;
        const double r  = std::sqrt (r2);
        const double rho = sigma0 * std::pow (r + rs, -0.5) * std::max (std::exp (-std::pow (r/r0, -xsi)), 1e-6);
        const double omegak2 = GM / std::pow (r + rs,  3.0) * NumHoles; // squared Keplerian for total M at center, is M total mass of binary?
        const double omega2  = omegak2;// * std::pow (1.0 + (3.0/16.0)*aBin*aBin/r2, 2.0); // missing pressure gradient term
              double vq = r * std::sqrt (omega2);
        const double pre = std::pow (vq / MachNumber, 2.0) * rho; // cs^2 * rho
        const double h0 = r / MachNumber;
        const double nu = ViscousAlpha * (vq / MachNumber) * h0; //ViscousAlpha * cs * h0
	    const double vr = -(3.0/2.0) * nu / r; //radial drift velocity

        vq -= r * omegaFrame; 

        const double vx = vq * (-y / r) + vr * (x / r);
        const double vy = vq * ( x / r) + vr * (y / r);

        return std::vector<double> {rho, vx, vy, 0, pre, x, y};
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
        const auto ag  = GravitationalAcceleration (x, y, t);

        auto S = StateArray();

        // Gravitational source term
        // ====================================================================
        S[DDD] = 0.0;
        S[S11] = dg * ag[0];
        S[S22] = dg * ag[1];
        S[S33] = 0.0;
        S[NRG] = 0.0; // dg * (ag[0] * vx + ag[1] * vy);
        S[XXX] = 0.0;
        S[YYY] = 0.0;

        // Fictitious forces due to rotating frame (without Euler force)

        S[S11] += dg * omegaFrame * ( 2.0 * vy + omegaFrame * x);
        S[S22] += dg * omegaFrame * (-2.0 * vx + omegaFrame * y);


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

        // Outer buffer
        // ====================================================================
        if (r > 0.5 * BufferRadius)
        {
            const double vk       = std::sqrt (-GravitationalPotential (x, y, t));
            const double torb     = 2.0 * M_PI * r / vk;
            const double tau      = torb * BetaBuffer / bufferZoneProfile (r);

            S[0] -= (state1.U[0] - state0.U[0]) / tau;
            S[1] -= (state1.U[1] - state0.U[1]) / tau;
            S[2] -= (state1.U[2] - state0.U[2]) / tau;
            S[3] -= (state1.U[3] - state0.U[3]) / tau;
            // S[4] -= (state1.U[4] - state0.U[4]) / tau;
        }

        // Sink
        // ====================================================================
        const double tsink = SinkTime (x, y, t);
        S[0] -= state1.U[0] / tsink;
        S[1] -= state1.U[1] / tsink;
        S[2] -= state1.U[2] / tsink;
        S[3] -= state1.U[3] / tsink;
        // S[4] -= state1.U[4] / tsink;

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

    if (MpiCommunicator::world().isThisMaster())
    {
        std::cout << user << std::endl;
    }

    SofteningRadius   = user["SofteningRadius"];
    BetaBuffer        = user["BetaBuffer"];
    BufferRadius      = user["BufferRadius"];
    ViscousAlpha      = user["ViscousAlpha"];
    MachNumber        = user["MachNumber"];
    NumHoles          = user["NumHoles"];
    aBin              = user["aBin"];
    SinkRadius        = user["SinkRadius"];
    InitialDataString = std::string (user["InitialData"]);

    if (InitialDataString == "LinearShear") initialData = initialDataLinearShear;
    if (InitialDataString == "Ring")        initialData = initialDataRing;
    if (InitialDataString == "Farris14")    initialData = initialDataFarris14;
    if (InitialDataString == "Tang17")      initialData = initialDataTang17;

    if (InitialDataString == "Tang17Rot")
    {       
        omegaFrame =  std::sqrt (2.0 * GM / (aBin * aBin * aBin));
        initialData = initialDataTang17Rot;
    }       
    auto logger    = std::make_shared<Logger>();
    auto writer    = std::make_shared<CheckpointWriter>();
    auto tseries   = std::make_shared<TimeSeriesManager>();
    auto scheduler = std::make_shared<TaskScheduler>();

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
    mg->setLowerUpper ({{-10.0, -10.0, 0.0}}, {{10.0, 10.0, 1.0}});
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

    cl->setGammaLawIndex (5. / 3);
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
    scheduler->schedule (taskRecomputeDt, TaskScheduler::Recurrence (0.0, 0.0, 1), "compute_dt"); // re-computing every time step

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
