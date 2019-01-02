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
static double BetaBuffer       = 1e-3;  // Orbital periods over which to relax to IC in outer buffer
static double ViscousAlpha     = 1e-1;  // Alpha viscosity parameter
static double VacuumCleaner    = 1.0;   // alpha_mini / alpha_disk (set to a number > 1 to simulate faster-than reasonable sinks)
static double MachNumber       = 10.0;  // 1/M = h/r
static int    NumHoles         = 2;     // Number of Black holes 1 or 2
static int    ScaleHeightFunc  = 0;     // 0 := H ~ cylindrical, 1 := H ~ min(r1,r2), 2 := H ~ interp(r1,r2)
static double GM               = 1.0;   // Newton G * mass of each component
static double aBin             = 1.0;   // Binary separation
static double SinkRadius       = 0.2;   // Sink radius
static double LdotEfficiency   = 1.0;   // Efficiency f to accrete L through the mini-disks (f = 1 is maximal)
static double DomainRadius     = 8.0;
static std::string InitialDataString = "Tang17";
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
        g1.x = 0.5 * aBin * std::cos (omegaBin * t);
        g1.y = 0.5 * aBin * std::sin (omegaBin * t);
        g2.x = 0.5 * aBin * std::cos (omegaBin * t + M_PI);
        g2.y = 0.5 * aBin * std::sin (omegaBin * t + M_PI);
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
        phi += -GM * std::pow (g.r * g.r + rs * rs, -0.5);
    }
    return phi;
}

static std::array<double, 2> GravitationalAcceleration (double x, double y, double t)
{
    const double rs = SofteningRadius;
    std::array<double, 2> a = {{0.0, 0.0}};

    for (const auto& g : SinkGeometries (x, y, t))
    {
        a[0] += -g.xhat * GM * g.r * std::pow (g.r * g.r + rs * rs, -1.5);
        a[1] += -g.yhat * GM * g.r * std::pow (g.r * g.r + rs * rs, -1.5);
    }
    return a;
}

/**
 * This function approximates the distance to either hole when that distance
 * is small, and the distance to the origin when it is large.
 */
double EffectiveRadius (double x, double y, double t)
{
    auto sinks = SinkGeometries (x, y, t);

    if (sinks.size() == 1)
    {
        return std::sqrt (x * x + y * y);
    }

    const double dx1 = x - sinks[0].x;
    const double dy1 = y - sinks[0].y;
    const double dx2 = x - sinks[1].x;
    const double dy2 = y - sinks[1].y;
    const double r1 = std::sqrt (dx1 * dx1 + dy1 * dy1);
    const double r2 = std::sqrt (dx2 * dx2 + dy2 * dy2);

    switch (ScaleHeightFunc)
    {
        case 0: return std::sqrt (x * x + y * y);
        case 1: return std::min (r1, r2);
        case 2: return (1.0 + 2.0 / (r1 / r2 + r2 / r1)) / (1.0 / r1 + 1.0 / r2);
    }
    throw;
}

static double SoundSpeedSquared (double x, double y, double t)
{
    const double phi = GravitationalPotential (x, y, t);
    const double cs2 = std::pow (MachNumber, -2) * -phi;
    return cs2;
}

static double SinkKernel (double r)
{
    const double omegaSink = std::sqrt (GM / std::pow (SinkRadius, 3));
    const double tvisc = 2. / 3 * MachNumber * MachNumber / ViscousAlpha / omegaSink;
    const double tsink = tvisc / VacuumCleaner;
    return r < SinkRadius ? tsink : HUGE_TIME_SCALE;
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
    // The tightness parameter below steepens the onset of the buffer. For a
    // value of 3 and when DomainRadius = 8, the return value is 1% at roughly
    // r = 7.2.
    const double tightness = 3;
    return 1.0 + std::tanh (tightness * (r - DomainRadius));
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
    double getSoundSpeedSquared (const double *P, double t) const;
};




// ============================================================================
ThinDiskNewtonianHydro::ThinDiskNewtonianHydro()
{
}

ConservationLaw::State ThinDiskNewtonianHydro::fromConserved (const Request& request, const double* U) const
{
    double P[7] = {0, 0, 0, 0, 0, 0, 0};

    P[RHO] = std::max (DENSITY_FLOOR, U[DDD]);
    P[PRE] = P[RHO] * getSoundSpeedSquared (U, request.simulationTime);
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
    const double cs = std::sqrt (getSoundSpeedSquared (P, request.simulationTime));
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
        case RHO: return "sigma";
        case V11: return "velocity1";
        case V22: return "velocity2";
        case V33: return "velocity3";
        case PRE: return "pressure";
        case XXX: return "x";
        case YYY: return "y";
        default: return "";
    }
}

double ThinDiskNewtonianHydro::getSoundSpeedSquared (const double* P, double t) const
{
    // Note: the argument is allowed to be P or U, since only x and y are read from it.
    const double x = P[XXX];
    const double y = P[YYY];
    return SoundSpeedSquared (x, y, t);
}

std::vector<double> ThinDiskNewtonianHydro::makeDiagnostics (const State& state) const
{
    const double t = status.simulationTime;
    const double x = state.P[XXX];
    const double y = state.P[YYY];
    const double dg = state.P[RHO];
    const double px = state.U[S11];
    const double py = state.U[S22];
    const double rs = SofteningRadius;
    const auto sinks = SinkGeometries (x, y, t);

    auto D = std::vector<double> (8);
    D[0] = dg;
    D[1] = x * py - y * px;

    if (sinks.size() >= 1)
    {
        const SinkGeometry g = sinks[0];
        const double agx = g.xhat * GM * g.r * std::pow (g.r * g.r + rs * rs, -1.5);
        const double agy = g.yhat * GM * g.r * std::pow (g.r * g.r + rs * rs, -1.5);
        const double tsink = SinkKernel (g.r);
        const double dgdot = dg / tsink;
        const double pxdot = px / tsink * LdotEfficiency;
        const double pydot = py / tsink * LdotEfficiency;
        const double Ldot1 = x * pydot - y * pxdot; // Total accretion torque from sink 1
        D[2] = dgdot;
        D[4] = dg * (g.x * agy - g.y * agx); // Gravitational torque on BH 1 (rb \times fg)
        D[6] = Ldot1;
    }
    if (sinks.size() >= 2)
    {
        const SinkGeometry g = sinks[1];
        const double agx = g.xhat * GM * g.r * std::pow (g.r * g.r + rs * rs, -1.5);
        const double agy = g.yhat * GM * g.r * std::pow (g.r * g.r + rs * rs, -1.5);
        const double tsink = SinkKernel (g.r);
        const double dgdot = dg / tsink;
        const double pxdot = px / tsink * LdotEfficiency;
        const double pydot = py / tsink * LdotEfficiency;
        const double Ldot2 = x * pydot - y * pxdot; // Total accretion torque from sink 2
        D[3] = dgdot;
        D[5] = dg * (g.x * agy - g.y * agx); // Gravitational torque on BH 2 (rb \times fg)
        D[7] = Ldot2;
    }
    return D;
}

std::vector<std::string> ThinDiskNewtonianHydro::getDiagnosticNames() const
{
    auto N = std::vector<std::string>(8);
    N[0] = "fluid_mass";
    N[1] = "fluid_angular_momentum";
    N[2] = "accretion_mdot_on_bh1";
    N[3] = "accretion_mdot_on_bh2";
    N[4] = "gravitational_torque_on_bh1";
    N[5] = "gravitational_torque_on_bh2";
    N[6] = "accretion_torque_on_bh1";
    N[7] = "accretion_torque_on_bh2";
    // N[7] = "accretion_torque_on_bh1-tang17";
    // N[8] = "accretion_torque_on_bh2-tang17";
    return N;
}




// ============================================================================
static void computeViscousFluxes2D (const Array& P, double dx, double dy, double t, Array& Fhat)
{
    for (int i = 0; i < P.shape()[0] - 1; ++i)
    {
        for (int j = 1; j < P.shape()[1] - 1; ++j)
        {
            const double dgL = P(i + 0, j, 0, RHO);
            const double dgR = P(i + 1, j, 0, RHO);
            const double xcL = P(i + 0, j, 0, XXX);
            const double xcR = P(i + 1, j, 0, XXX);
            const double ycL = P(i + 0, j, 0, YYY);
            const double ycR = P(i + 1, j, 0, YYY);
            const double uxL0 = P(i + 0, j - 1, 0, V11);
            const double uxL1 = P(i + 0, j + 0, 0, V11);
            const double uxL2 = P(i + 0, j + 1, 0, V11);
            const double uxR0 = P(i + 1, j - 1, 0, V11);
            const double uxR1 = P(i + 1, j + 0, 0, V11);
            const double uxR2 = P(i + 1, j + 1, 0, V11);
            const double uyL0 = P(i + 0, j - 1, 0, V22);
            const double uyL1 = P(i + 0, j + 0, 0, V22);
            const double uyL2 = P(i + 0, j + 1, 0, V22);
            const double uyR0 = P(i + 1, j - 1, 0, V22);
            const double uyR1 = P(i + 1, j + 0, 0, V22);
            const double uyR2 = P(i + 1, j + 1, 0, V22);
            const double dg = 0.5 * (dgL + dgR);
            const double xc = 0.5 * (xcL + xcR);
            const double yc = 0.5 * (ycL + ycR);
            const double cs = std::sqrt (SoundSpeedSquared (xc, yc, t));
            const double rt = EffectiveRadius (xc, yc, t);
            const double nu = ViscousAlpha * cs * rt / MachNumber;
            const double mu = dg * nu;

            const double dxux = (uxR1 - uxL1) / dx;
            const double dxuy = (uyR1 - uyL1) / dx;
            const double dyux = (uxR2 - uxR0 + uxL2 - uxL0) / (4 * dy);
            const double dyuy = (uyR2 - uyR0 + uyL2 - uyL0) / (4 * dy);
            const double tauxx = mu * (dxux - dyuy);
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
            const double xcL = P(i, j + 0, 0, XXX);
            const double xcR = P(i, j + 1, 0, XXX);
            const double ycL = P(i, j + 0, 0, YYY);
            const double ycR = P(i, j + 1, 0, YYY);
            const double uxL0 = P(i - 1, j + 0, 0, V11);
            const double uxL1 = P(i + 0, j + 0, 0, V11);
            const double uxL2 = P(i + 1, j + 0, 0, V11);
            const double uxR0 = P(i - 1, j + 1, 0, V11);
            const double uxR1 = P(i + 0, j + 1, 0, V11);
            const double uxR2 = P(i + 1, j + 1, 0, V11);
            const double uyL0 = P(i - 1, j + 0, 0, V22);
            const double uyL1 = P(i + 0, j + 0, 0, V22);
            const double uyL2 = P(i + 1, j + 0, 0, V22);
            const double uyR0 = P(i - 1, j + 1, 0, V22);
            const double uyR1 = P(i + 0, j + 1, 0, V22);
            const double uyR2 = P(i + 1, j + 1, 0, V22);
            const double dg = 0.5 * (dgL + dgR);
            const double xc = 0.5 * (xcL + xcR);
            const double yc = 0.5 * (ycL + ycR);
            const double cs = std::sqrt (SoundSpeedSquared (xc, yc, t));
            const double rt = EffectiveRadius (xc, yc, t);
            const double nu = ViscousAlpha * cs * rt / MachNumber;
            const double mu = dg * nu;

            const double dyux = (uxR1 - uxL1) / dy;
            const double dyuy = (uyR1 - uyL1) / dy;
            const double dxux = (uxR2 - uxR0 + uxL2 - uxL0) / (4 * dx);
            const double dxuy = (uyR2 - uyR0 + uyL2 - uyL0) / (4 * dx);
            const double tauyx = mu * (dyux + dxuy);
            const double tauyy = mu * (dyuy - dxux);

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
    return [dx,dy] (const Cow::Array& P, Cow::Array& F, double t)
    {
        computeViscousFluxes2D (P, dx, dy, t, F);
    };
}




// ============================================================================
int BinaryTorque::run (int argc, const char* argv[])
{
    auto user = Variant::NamedValues();

    user["outdir"]  = "data";
    user["restart"] = "";
    user["tfinal"]  = 16.0;
    user["serial"]  = 0;
    user["cpi"]     = 0.25;
    user["cpf"]     = "single"; // or multiple
    user["tsi"]     = 0.1;
    user["cfl"]     = 0.3;
    user["plm"]     = 1.5;
    user["N"]       = 128;


    user["SofteningRadius"]  = SofteningRadius;
    user["BetaBuffer"]       = BetaBuffer;
    user["DomainRadius"]     = DomainRadius;
    user["ViscousAlpha"]     = ViscousAlpha;
    user["MachNumber"]       = MachNumber;
    user["NumHoles"]         = NumHoles;
    user["ScaleHeightFunc"]  = ScaleHeightFunc;
    user["aBin"]             = aBin;
    user["SinkRadius"]       = SinkRadius;
    user["LdotEfficiency"]   = LdotEfficiency;
    user["InitialData"]      = InitialDataString;
    user["VacuumCleaner"]    = VacuumCleaner;


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

    auto initialDataZrake = [&] (double x, double y, double z) -> std::vector<double>
    {
        const auto   phi        = GravitationalPotential (x, y, status.simulationTime);
        const double cs2        = SoundSpeedSquared      (x, y, status.simulationTime);
        const double sigma0     = 1.0;
        const double r0         = 2.5 * aBin;
        const double rs         = SofteningRadius;
        const double r2         = x * x + y * y;
        const double r          = std::sqrt (r2);
        const double sigma      = std::pow (r/r0 + rs/r0, -0.5) * sigma0;
        const double vq         = std::sqrt (-phi);
        const double pre        = cs2 * sigma;
        const double vx         = vq * (-y / r);
        const double vy         = vq * ( x / r);
        return std::vector<double> {sigma, vx, vy, 0, pre, x, y};
    };

    auto initialDataTang17 = [&] (double x, double y, double z) -> std::vector<double>
    {
        // Initial conditions from Tang+ (2017) MNRAS 469, 4258
        // MISSING?: radial drift velocity, pressure gradient for omegak2
        // Check Yike/Chris Disco setup
        const auto   ag         = GravitationalAcceleration (x, y, status.simulationTime);
        const double rs         = SofteningRadius;
        const double r0         = 2.5 * aBin;
        const double xsi        = 10.0;
        const double sigma0     = 1.0;
        const double r2         = x * x + y * y;
        const double r          = std::sqrt (r2);
        const double cs2        = SoundSpeedSquared (x, y, status.simulationTime);
        const double cs2Deriv   = +(ag[0] * x + ag[1] * y) / r * std::pow (MachNumber, -2);
        const double sigma      = sigma0 * std::pow (r/r0 + rs/r0, -0.5) * std::max (std::exp (-std::pow (r / r0, -xsi)), 1e-6);
        const double sigmaDeriv = sigma0 * std::pow (r/r0 + rs/r0, -1.5) * -0.5 / r0; // neglects the cavity cutoff ^^
        const double dPdr       = cs2 * sigmaDeriv + cs2Deriv * sigma;
        const double omega2     = r < r0 ? GM / (4 * r0) : -(ag[0] * x + ag[1] * y) / r2 + dPdr / (sigma * r); // Not sure why GM / (4 * r0)... check this
        const double vq         = r * std::sqrt (omega2); // defined to be in the inertial frame
        const double pre        = cs2 * sigma;
        const double h0         = r / MachNumber;
        const double nu         = ViscousAlpha * std::sqrt (cs2) * h0; // ViscousAlpha * cs * h0
        const double vr         = -(3.0 / 2.0) * nu / (r + rs); // radial drift velocity (CHECK)
        const double vx         = vq * (-y / r) + vr * (x / r);
        const double vy         = vq * ( x / r) + vr * (y / r);

        return std::vector<double> {sigma, vx, vy, 0, pre, x, y};
    };

    // Source terms
    // ------------------------------------------------------------------------
    auto sourceTermsFunction = [&] (double x, double y, double z, double t, StateArray P)
    {
        const double dg = P[0];
        const double r  = std::sqrt (x * x + y * y);
        const auto ag   = GravitationalAcceleration (x, y, t);
        auto S = StateArray();


        // Gravitational source term
        // ====================================================================
        S[DDD] = 0.0;
        S[S11] = dg * ag[0];
        S[S22] = dg * ag[1];
        S[S33] = 0.0;
        S[NRG] = 0.0;
        S[XXX] = 0.0;
        S[YYY] = 0.0;

        ConservationLaw::Request request;
        request.getPrimitive = true;
        request.getConserved = true;
        request.getFluxes = false;
        request.getEigenvalues = false;
        request.areaElement[0] = 1.0;
        request.areaElement[1] = 0.0;
        request.areaElement[2] = 0.0;
        request.simulationTime = t;

        auto state0 = cl->fromPrimitive (request, &initialData (x, y, z)[0]);
        auto state1 = cl->fromPrimitive (request, &P[0]);


        // Outer buffer
        // ====================================================================
        const double vk       = std::sqrt (-GravitationalPotential (x, y, t));
        const double torb     = 2.0 * M_PI * r / vk;
        const double tau      = torb * BetaBuffer / bufferZoneProfile (r);

        S[0] -= (state1.U[0] - state0.U[0]) / tau;
        S[1] -= (state1.U[1] - state0.U[1]) / tau;
        S[2] -= (state1.U[2] - state0.U[2]) / tau;


        // Sink
        // ====================================================================
        const double tsink = SinkTime (x, y, t);
        S[0] -= state1.U[0] / tsink;
        S[1] -= state1.U[1] / tsink * LdotEfficiency;
        S[2] -= state1.U[2] / tsink * LdotEfficiency;

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
    DomainRadius      = user["DomainRadius"];
    ViscousAlpha      = user["ViscousAlpha"];
    MachNumber        = user["MachNumber"];
    NumHoles          = user["NumHoles"];
    ScaleHeightFunc   = user["ScaleHeightFunc"];
    aBin              = user["aBin"];
    SinkRadius        = user["SinkRadius"];
    LdotEfficiency    = user["LdotEfficiency"];
    VacuumCleaner     = user["VacuumCleaner"];
    InitialDataString = std::string (user["InitialData"]);


    if (InitialDataString == "LinearShear") initialData = initialDataLinearShear;
    if (InitialDataString == "Ring")        initialData = initialDataRing;
    if (InitialDataString == "Zrake")       initialData = initialDataZrake;
    if (InitialDataString == "Tang17")      initialData = initialDataTang17;


    auto logger    = std::make_shared<Logger>();
    auto writer    = std::make_shared<CheckpointWriter>();
    auto tseries   = std::make_shared<TimeSeriesManager>();
    auto scheduler = std::make_shared<TaskScheduler>();

    auto bd = std::shared_ptr<BlockDecomposition>();
    auto mg = std::shared_ptr<MeshGeometry> (new CartesianMeshGeometry);
    auto rs = std::shared_ptr<RiemannSolver> (new HllcNewtonianHydroRiemannSolver);
    auto fs = std::shared_ptr<IntercellFluxScheme> (new MethodOfLinesPlm);
    auto bc = std::shared_ptr<BoundaryCondition>(new OutflowBoundaryCondition);


    // Set up grid shape.
    // ------------------------------------------------------------------------
    int N = int (user["N"]);
    auto cs = Shape {{ N, N, 1, 1, 1 }};
    auto bs = Shape {{ 2, 2, 0, 0, 0 }};

    mg->setCellsShape (cs);
    mg->setLowerUpper ({{-DomainRadius, -DomainRadius, 0.0}}, {{DomainRadius, DomainRadius, 1.0}});

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
    auto advance   = [&] (double dt) { return ss->advance (*md, status.simulationTime, dt); };
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

    auto taskTimeSeries = [&] (SimulationStatus, int rep)
    {
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

    scheduler->schedule (taskCheckpoint, TaskScheduler::Recurrence (user["cpi"]), "checkpoint");
    scheduler->schedule (taskRecomputeDt, TaskScheduler::Recurrence (0.0, 0.0, 1), "compute_dt");
    scheduler->schedule (taskTimeSeries, TaskScheduler::Recurrence (0.0, 0.0, 1), "time_series");

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
