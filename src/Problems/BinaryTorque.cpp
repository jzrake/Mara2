#include <mpi.h>
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
static double SofteningRadius  =  0.2;  // rs^2, where Fg = G M m / (r^2 + rs^2). If rs < 0, then rs is measured in grid cells
static double BetaBuffer       = 1e-3;  // Orbital periods over which to relax to IC in outer buffer
static double ViscousAlpha     = 1e-1;  // Alpha viscosity parameter
static double VacuumCleaner    = 1.0;   // if >0 =alpha_mini/alpha_disk (> 1 faster-than-viscous sinks) if<0 tsink=tBinary/abs(VacuuumCleaer)
static double MachNumber       = 10.0;  // 1/M = h/r
static int    NumHoles         = 2;     // Number of Black holes 1 or 2
static int    ScaleHeightFunc  = 0;     // 0 := H ~ cylindrical, 1 := H ~ min(r1,r2), 2 := H ~ interp(r1,r2)
static int    SinkKernelFunc   = 1;     // 0 := tophat, 1 := Gaussian 2 := super-Gaussian
static int    LiveBinary       = 0;     // Move the binary! (enables RK3)
static int    LocalViscousSink = 0;     // Use local viscous time scale for sink
static double GM               = 1.0;   // Newton G * mass of binary
static double aBin             = 1.0;   // Binary separation
static double BinaryMassRatio  = 1.0;   // Binary mass ratio q = m2/m1
static double DiskMassRatio    = 0.1;   // sigma0 = DiskMassRatio * M / aBin^2, where sigma0 is used in Tang17
static double Eccentricity     = 0.0;   // Binary eccentricity e
static double SinkRadius       = 0.2;   // Sink radius
static double LdotEfficiency   = 1.0;   // Efficiency f to accrete L through the mini-disks (f = 1 is maximal)
static double DomainRadius     = 8.0;
static std::string VersionNumber     = "unknown";
static std::string InitialDataString = "Tang17";
static SimulationStatus status;
static std::array<double, 8> globalStarParticleData;
static double vStarMax = 0.0; // is written to in the sinkParticleLocations function




// ============================================================================
struct SinkGeometry
{
    double x;
    double y;
    double r;
    double xhat;
    double yhat;
};


static std::array<double, 4> SinkLocations (std::array<double, 8> T)
{
    return {T[0], T[1], T[2], T[3]};
}

static std::vector<SinkGeometry> SinkGeometries (double x, double y, std::array<double, 8> T)
{
    if (NumHoles == 1)
    {
        SinkGeometry g1, g2;

        g1.x = 0.0;
        g1.y = 0.0;
        g2.x = 0.0;
        g2.y = 0.0;

        g1.r = std::sqrt ((x - g1.x) * (x - g1.x) + (y - g1.y) * (y - g1.y));
        g2.r = std::sqrt ((x - g2.x) * (x - g2.x) + (y - g2.y) * (y - g2.y));
        g1.xhat = (x - g1.x) / g1.r;
        g1.yhat = (y - g1.y) / g1.r;
        g2.xhat = (x - g2.x) / g2.r;
        g2.yhat = (y - g2.y) / g2.r;

        return {g1, g2};
    }
    if (NumHoles == 2)
    {
        const auto sinkLocations = SinkLocations (T);

        SinkGeometry g1, g2;

        g1.x = sinkLocations[0];
        g1.y = sinkLocations[1];
        g2.x = sinkLocations[2];
        g2.y = sinkLocations[3];

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

static double GravitationalPotential (double x, double y, std::array<double, 8> T)
{
    const double rs = SofteningRadius;
    const double GM1 = GM  / (1.0 + BinaryMassRatio);
    const double GM2 = GM1 * BinaryMassRatio;

    auto s = SinkGeometries (x, y, T);

    double phi  = -GM1 * std::pow (s[0].r * s[0].r + rs * rs, -0.5);
           phi += -GM2 * std::pow (s[1].r * s[1].r + rs * rs, -0.5);
    return phi;
}

static std::array<double, 2> GravitationalAcceleration (double x, double y, std::array<double, 8> T)
{
    std::array<double, 2> a = {{0.0, 0.0}};
    const double rs = SofteningRadius;
    const double GM1 = GM / (1.0 + BinaryMassRatio);
    const double GM2 = GM1 *  BinaryMassRatio;

    auto s = SinkGeometries (x, y, T);

    a[0] += -s[0].xhat * GM1 * s[0].r * std::pow (s[0].r * s[0].r + rs * rs, -1.5);
    a[0] += -s[1].xhat * GM2 * s[1].r * std::pow (s[1].r * s[1].r + rs * rs, -1.5);

    a[1] += -s[0].yhat * GM1 * s[0].r * std::pow (s[0].r * s[0].r + rs * rs, -1.5);
    a[1] += -s[1].yhat * GM2 * s[1].r * std::pow (s[1].r * s[1].r + rs * rs, -1.5);

    return a;
}

/**
 * This function approximates the distance to either hole when that distance
 * is small, and the distance to the origin when it is large.
 */
double EffectiveRadius (double x, double y, std::array<double, 8> T)
{
    auto sinks = SinkGeometries (x, y, T);

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
    const double r  = std::sqrt (x * x + y * y);

    switch (ScaleHeightFunc)
    {
        case 0: return r;
        case 1: return std::min (r1, r2);
        case 2: return (1.0 + 2.0 / (r1 / r2 + r2 / r1)) / (1.0 / r1 + 1.0 / r2);
        case 3: return r < 2.5 * aBin ? 0.0 : r;
    }
    throw;
}

static double SoundSpeedSquared (double x, double y, std::array<double, 8> T)
{
    const double phi = GravitationalPotential (x, y, T);
    const double cs2 = std::pow (MachNumber, -2) * -phi;
    return cs2;
}

static double SinkKernel1 (double r)
{
    const double sr2 = SinkRadius * SinkRadius;
    double tsink;

    if (VacuumCleaner > 0)
    {
        const double GM1 = GM / (1.0 + BinaryMassRatio);
        const double rs = (LocalViscousSink > 0) ? r : SinkRadius;
        const double omegaSink = std::sqrt (GM1 / std::pow (rs, 3));
        const double tvisc = 2.0 / 3.0 * MachNumber * MachNumber / ViscousAlpha / omegaSink;
        tsink = tvisc / VacuumCleaner;
    }
    else
    {
        const double omegaBin = aBin == 0.0 ? 0.0 : std::sqrt (GM / (aBin * aBin * aBin));
        const double tBinary = 2.0 * M_PI / omegaBin;
        tsink = tBinary / std::abs(VacuumCleaner);
    }
    switch (SinkKernelFunc)
    {
        case 0: return r < SinkRadius ? 1.0 / tsink : 0.0;
        case 1: return (1.0 / tsink) * std::exp(-0.5 * r * r / sr2);
        case 2: return (1.0 / tsink) * std::exp(-0.5 * r * r * r * r / sr2 / sr2);
        default: throw std::invalid_argument("SinkKernelFunc must be 0, 1, or 2");
    }
}

static double SinkKernel2 (double r)
{
    const double sr2 = SinkRadius * SinkRadius;
    double tsink;

    if (VacuumCleaner > 0)
    {
        const double GM2 = GM * BinaryMassRatio / (1.0 + BinaryMassRatio);
        const double rs = (LocalViscousSink > 0) ? r : SinkRadius;
        const double omegaSink = std::sqrt (GM2 / std::pow (rs, 3));
        const double tvisc = 2.0 / 3.0 * MachNumber * MachNumber / ViscousAlpha / omegaSink;
        tsink = tvisc / VacuumCleaner;
    }
    else
    {
        const double omegaBin = aBin == 0.0 ? 0.0 : std::sqrt (GM / (aBin * aBin * aBin));
        const double tBinary = 2.0 * M_PI / omegaBin;
        tsink = tBinary / std::abs(VacuumCleaner);
    }
    switch (SinkKernelFunc)
    {
        case 0: return r < SinkRadius ? 1.0 / tsink : 0.0;
        case 1: return (1.0 / tsink) * std::exp(-0.5 * r * r / sr2);
        case 2: return (1.0 / tsink) * std::exp(-0.5 * r * r * r * r / sr2 / sr2);
        default: throw std::invalid_argument("SinkKernelFunc must be 0, 1, or 2");
    }
}

static double SinkBeta (double x, double y, std::array<double, 8> T)
{
    auto s = SinkGeometries (x, y, T);
    double beta2 = 0.0;

    double const beta1 = SinkKernel1 (s[0].r);

    if (s.size() > 1)
    {
        beta2 = SinkKernel2 (s[1].r);
    }
    return std::max(beta1, beta2);
}


// static double SinkTime (double x, double y, std::array<double, 8> T)
// {
//     double tmin = HUGE_TIME_SCALE;
//     auto s = SinkGeometries (x, y, T);

//     if (SinkKernel1 (s[0].r) < tmin)
//     {
//         tmin = SinkKernel1 (s[0].r);
//     }

//     if (SinkKernel2 (s[1].r) < tmin)
//     {
//         tmin = SinkKernel2 (s[1].r);
//     }

//     return tmin;
// }

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
    double getSoundSpeedSquared (const double *P, std::array<double, 8> T) const;
};

// ============================================================================
ThinDiskNewtonianHydro::ThinDiskNewtonianHydro()
{
}

ConservationLaw::State ThinDiskNewtonianHydro::fromConserved (const Request& request, const double* U) const
{
    double P[7] = {0, 0, 0, 0, 0, 0, 0};

    P[RHO] = std::max (DENSITY_FLOOR, U[DDD]);
    P[PRE] = P[RHO] * getSoundSpeedSquared (U, request.auxiliaryData);
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
    const double cs = std::sqrt (getSoundSpeedSquared (P, request.auxiliaryData));
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

double ThinDiskNewtonianHydro::getSoundSpeedSquared (const double* P, std::array<double, 8> T) const
{
    // Note: the argument is allowed to be P or U, since only x and y are read from it.
    const double x = P[XXX];
    const double y = P[YYY];
    return SoundSpeedSquared (x, y, T);
}

std::vector<double> ThinDiskNewtonianHydro::makeDiagnostics (const State& state) const
{
    const auto T = globalStarParticleData;
    const double x = state.P[XXX];
    const double y = state.P[YYY];
    const double dg = state.P[RHO];
    const double px = state.U[S11];
    const double py = state.U[S22];
    const double rs = SofteningRadius;
    const double GM1 = GM / (1.0 + BinaryMassRatio);
    const double GM2 = GM1 *  BinaryMassRatio;
    const auto sinks = SinkGeometries (x, y, T);

    auto D = std::vector<double> (18);
    D[0] = dg;
    D[1] = x * py - y * px;

    if (sinks.size() >= 1)
    {
        const SinkGeometry g = sinks[0];
        const double agx = g.xhat * GM1 * g.r * std::pow (g.r * g.r + rs * rs, -1.5);
        const double agy = g.yhat * GM1 * g.r * std::pow (g.r * g.r + rs * rs, -1.5);
        const double betasink = SinkKernel1 (g.r);
        const double dgdot = dg * betasink;
        const double pxdot = px * betasink * LdotEfficiency;
        const double pydot = py * betasink * LdotEfficiency;
        const double Ldot1 = x * pydot - y * pxdot; // Total accretion torque from sink 1
        const double deltax = x - g.x;
        const double deltay = y - g.y;
        const double sdot1 = pydot * deltax - pxdot * deltay;
        D[2] = dgdot;
        D[4] = dg * (g.x * agy - g.y * agx); // Gravitational torque on BH 1 (rb \times fg)
        D[6] = Ldot1;
        D[8] = sdot1;
    }
    if (sinks.size() >= 2)
    {
        const SinkGeometry g = sinks[1];
        const double agx = g.xhat * GM2 * g.r * std::pow (g.r * g.r + rs * rs, -1.5);
        const double agy = g.yhat * GM2 * g.r * std::pow (g.r * g.r + rs * rs, -1.5);
        const double betasink = SinkKernel2 (g.r);
        const double dgdot = dg * betasink;
        const double pxdot = px * betasink * LdotEfficiency;
        const double pydot = py * betasink * LdotEfficiency;
        const double Ldot2 = x * pydot - y * pxdot; // Total accretion torque from sink 2
        const double deltax = x - g.x;
        const double deltay = y - g.y;
        const double sdot2 = pydot * deltax - pxdot * deltay;
        D[3] = dgdot;
        D[5] = dg * (g.x * agy - g.y * agx); // Gravitational torque on BH 2 (rb \times fg)
        D[7] = Ldot2;
        D[9] = sdot2;
    }

    D[10] = globalStarParticleData[0];
    D[11] = globalStarParticleData[1];
    D[12] = globalStarParticleData[2];
    D[13] = globalStarParticleData[3];
    D[14] = globalStarParticleData[4];
    D[15] = globalStarParticleData[5];
    D[16] = globalStarParticleData[6];
    D[17] = globalStarParticleData[7];

    return D;
}

std::vector<std::string> ThinDiskNewtonianHydro::getDiagnosticNames() const
{
    auto N = std::vector<std::string>(18);
    N[0] = "fluid_mass";
    N[1] = "fluid_angular_momentum";
    N[2] = "accretion_mdot_on_bh1";
    N[3] = "accretion_mdot_on_bh2";
    N[4] = "gravitational_torque_on_bh1";
    N[5] = "gravitational_torque_on_bh2";
    N[6] = "accretion_torque_on_bh1";
    N[7] = "accretion_torque_on_bh2";
    N[8] = "accretion_sdot_on_bh1";
    N[9] = "accretion_sdot_on_bh2";
    N[10] = "bh1_x";
    N[11] = "bh1_y";
    N[12] = "bh2_x";
    N[13] = "bh2_y";
    N[14] = "bh1_vx";
    N[15] = "bh1_vy";
    N[16] = "bh2_vx";
    N[17] = "bh2_vy";
    return N;
}




// ============================================================================
static void computeViscousFluxes2D (const Array& P, double dx, double dy, std::array<double, 8> T, Array& Fhat)
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
            const double cs = std::sqrt (SoundSpeedSquared (xc, yc, T));
            const double rt = EffectiveRadius (xc, yc, T);
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
            const double cs = std::sqrt (SoundSpeedSquared (xc, yc, T));
            const double rt = EffectiveRadius (xc, yc, T);
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
    return [dx,dy] (const Cow::Array& P, Cow::Array& F, std::array<double, 8> t)
    {
        computeViscousFluxes2D (P, dx, dy, t, F);
    };
}

static std::vector<double> initialStarParticleLocations()
{
    const double omegaBin = aBin == 0.0 ? 0.0 : std::sqrt (GM / (aBin * aBin * aBin));

    return {
        +0.5 * aBin, 0.0,
        -0.5 * aBin, 0.0,
        0.0, +omegaBin * 0.5 * aBin,
        0.0, -omegaBin * 0.5 * aBin,
    };
};

static std::vector<double> starParticleLocations (double t)
{
    static double Tolerance = 1e-6;

    std::vector<double> sinkLocations =
    {
        0, 0, 0, 0, // positions
        0, 0, 0, 0, // velocities
    };

    const double e = Eccentricity;
    const double q = BinaryMassRatio;
    const double omegaBin = aBin == 0.0 ? 0.0 : std::sqrt (GM / (aBin * aBin * aBin));
    const double mu = q / (1.0 + q);
    double v1, v2;

    if (Eccentricity > 0.0)
    {
        double M = omegaBin * t; //Mean anomoly
        // minimize f(E) = E - e * sinE - M using Newton-Rapheson
        double    E = M; // eccentric anomoly - set first guess to M
        double fofE = E - e * std::sin(E) - M;

        while (std::abs(fofE) > Tolerance)
        {
            double dfdE = 1.0 - e * std::cos(E);
            E   -= fofE/dfdE;
            fofE = E - e * std::sin(E) - M;
        }
        double x = aBin * (std::cos(E) - e);
        double y = aBin * std::sqrt(1.0 - e * e) * std::sin(E);
        double r = std::sqrt(x*x + y*y);
        double f = std::atan2(y,x);

        sinkLocations[0] = mu *     r * std::cos(f);
        sinkLocations[1] = mu *     r * std::sin(f);
        sinkLocations[2] = mu / q * r * std::cos(f + M_PI);
        sinkLocations[3] = mu / q * r * std::sin(f + M_PI);

        double rdot  = omegaBin * aBin / std::sqrt(1 - e * e) * e * std::sin(f);
        double rfdot = omegaBin * aBin / std::sqrt(1 - e * e) * (1.0 + e * std::cos(f));

        v1 = mu     * std::sqrt(rdot * rdot + rfdot * rfdot);
        v2 = mu / q * std::sqrt(rdot * rdot + rfdot * rfdot);
    }
    else
    {
        sinkLocations[0] = mu     * aBin * std::cos (omegaBin * t);
        sinkLocations[1] = mu     * aBin * std::sin (omegaBin * t);
        sinkLocations[2] = mu / q * aBin * std::cos (omegaBin * t + M_PI);
        sinkLocations[3] = mu / q * aBin * std::sin (omegaBin * t + M_PI);
        v1 = mu     * aBin * omegaBin;
        v2 = mu / q * aBin * omegaBin;
    }

    vStarMax = std::max (v1, v2);

    return sinkLocations;
}

static std::vector<double> starParticleDerivatives (const Cow::Array& P, const std::vector<double>& starParticles)
{
    double GM1 = GM  / (1.0 + BinaryMassRatio);
    double GM2 = GM1 * BinaryMassRatio;
    const double rs = SofteningRadius;

    const double x1 = starParticles[0];
    const double y1 = starParticles[1];
    const double x2 = starParticles[2];
    const double y2 = starParticles[3];
    const double vx1 = starParticles[4];
    const double vy1 = starParticles[5];
    const double vx2 = starParticles[6];
    const double vy2 = starParticles[7];
    const double x1dot = vx1;
    const double y1dot = vy1;
    const double x2dot = vx2;
    const double y2dot = vy2;

    const double x12 = x2 - x1;
    const double y12 = y2 - y1;
    const double r2 = x12 * x12 + y12 * y12;
    const double x12hat = x12 / std::sqrt (r2);
    const double y12hat = y12 / std::sqrt (r2);

    const double vx1dot =  GM2 / r2 * x12hat;
    const double vy1dot =  GM2 / r2 * y12hat;
    const double vx2dot = -GM1 / r2 * x12hat;
    const double vy2dot = -GM1 / r2 * y12hat;

    // WARNING: assuming two guard zones here:
    const double dx = P(3, 0, 0, XXX) - P(2, 0, 0, XXX);
    const double dy = P(0, 3, 0, YYY) - P(0, 2, 0, YYY);

    double vx1dotGas = 0.0;
    double vy1dotGas = 0.0;
    double vx2dotGas = 0.0;
    double vy2dotGas = 0.0;

    for (int i = 2; i < P.shape()[0] - 2; ++i)
    {
        for (int j = 2; j < P.shape()[1] - 2; ++j)
        {
            const double sigma = P(i, j, 0, RHO);
            const double mc = sigma * dx * dy;
            const double xc = P(i, j, 0, XXX);
            const double yc = P(i, j, 0, YYY);
            const double x1c = xc - x1;
            const double y1c = yc - y1;
            const double x2c = xc - x2;
            const double y2c = yc - y2;
            const double r1csqu = x1c * x1c + y1c * y1c;
            const double r2csqu = x2c * x2c + y2c * y2c;

            vx1dotGas += mc * std::pow (r1csqu + rs * rs, -1.5) * x1c;
            vy1dotGas += mc * std::pow (r1csqu + rs * rs, -1.5) * y1c;
            vx2dotGas += mc * std::pow (r2csqu + rs * rs, -1.5) * x2c;
            vy2dotGas += mc * std::pow (r2csqu + rs * rs, -1.5) * y2c;
        }
    }

    // agas is the gravitational acceleration on the star from the gas.
    std::vector<double> agasLocal = {vx1dotGas, vy1dotGas, vx2dotGas, vy2dotGas};
    std::vector<double> agasTotal = MpiCommunicator::world().sum (agasLocal);

    return {
        x1dot, y1dot,
        x2dot, y2dot,
        vx1dot + agasTotal[0], vy1dot + agasTotal[1],
        vx2dot + agasTotal[2], vy2dot + agasTotal[3],
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
    user["tspurge"] = 0; // set to 1 if you want the time series data to be purged after each checkpoint

    user["SofteningRadius"]  = SofteningRadius;
    user["BetaBuffer"]       = BetaBuffer;
    user["DomainRadius"]     = DomainRadius;
    user["ViscousAlpha"]     = ViscousAlpha;
    user["MachNumber"]       = MachNumber;
    user["NumHoles"]         = NumHoles;
    user["ScaleHeightFunc"]  = ScaleHeightFunc;
    user["SinkKernelFunc"]   = SinkKernelFunc;
    user["aBin"]             = aBin;
    user["BinaryMassRatio"]  = BinaryMassRatio;
    user["DiskMassRatio"]    = DiskMassRatio;
    user["Eccentricity"]     = Eccentricity;
    user["SinkRadius"]       = SinkRadius;
    user["LocalViscousSink"] = LocalViscousSink;
    user["LiveBinary"]       = LiveBinary;
    user["LdotEfficiency"]   = LdotEfficiency;
    user["InitialData"]      = InitialDataString;
    user["VacuumCleaner"]    = VacuumCleaner;
    user["VersionNumber"]    = VersionNumber;


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
        const double r2 = x * x + y * y;
        const double r  = std::sqrt (r2);
        const double rho = 0.1 + std::exp (-std::pow (r - 4, 2));
        const double vq = std::sqrt (-GravitationalPotential (x, y, globalStarParticleData));
        const double pre = std::pow (vq / MachNumber, 2.0) * rho; // cs^2 * rho
        const double qhX = -y / r;
        const double qhY =  x / r;
        return std::vector<double> {rho, vq * qhX, vq * qhY, 0, pre, x, y};
    };

    // auto initialDataZrake = [&] (double x, double y, double z) -> std::vector<double>
    // {
    //     const auto   phi        = GravitationalPotential (x, y, globalStarParticleData);
    //     const double cs2        = SoundSpeedSquared      (x, y, globalStarParticleData);
    //     const double sigma0     = 1.0;
    //     const double r0         = 2.5 * aBin;
    //     const double rs         = SofteningRadius;
    //     const double r2         = x * x + y * y;
    //     const double r          = std::sqrt (r2);
    //     const double sigma      = std::pow (r/r0 + rs/r0, -0.5) * sigma0;
    //     const double vq         = std::sqrt (-phi);
    //     const double pre        = cs2 * sigma;
    //     const double vx         = vq * (-y / r);
    //     const double vy         = vq * ( x / r);
    //     return std::vector<double> {sigma, vx, vy, 0, pre, x, y};
    // };

    auto initialDataTang17 = [&] (double x, double y, double z) -> std::vector<double>
    {
        // Initial conditions from Tang+ (2017) MNRAS 469, 4258
        // MISSING?: radial drift velocity, pressure gradient for omegak2
        // Check Yike/Chris Disco setup
        const auto   T          = globalStarParticleData;
        const auto   ag         = GravitationalAcceleration (x, y, T);
        const double rs         = SofteningRadius;
        const double r0         = 2.5 * aBin;
        const double xsi        = 10.0;
        const double sigma0     = DiskMassRatio * GM / aBin / aBin; // NOTE: there's a G here, meaning G = 1 throughout
        const double r2         = x * x + y * y;
        const double r          = std::sqrt (r2);
        const double cs2        = SoundSpeedSquared (x, y, T);
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
    auto sourceTermsFunction = [&] (double x, double y, double z, std::array<double, 8> t, StateArray P)
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
        request.auxiliaryData = t;

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
        const double betasink = SinkBeta (x, y, t);
        S[0] -= state1.U[0] * betasink;
        S[1] -= state1.U[1] * betasink * LdotEfficiency;
        S[2] -= state1.U[2] * betasink * LdotEfficiency;

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
    SinkKernelFunc    = user["SinkKernelFunc"];
    aBin              = user["aBin"];
    BinaryMassRatio   = user["BinaryMassRatio"];
    DiskMassRatio     = user["DiskMassRatio"];
    Eccentricity      = user["Eccentricity"];
    SinkRadius        = user["SinkRadius"];
    LocalViscousSink  = user["LocalViscousSink"];
    LiveBinary        = user["LiveBinary"];
    LdotEfficiency    = user["LdotEfficiency"];
    VacuumCleaner     = user["VacuumCleaner"];
    InitialDataString = std::string (user["InitialData"]);
    VersionNumber     = std::string (user["VersionNumber"]);

    if (SofteningRadius < 0.0)
    {
        double ncells = -SofteningRadius;
        SofteningRadius = ncells * 2.0 * DomainRadius / int (user["N"]);
    }

    if (InitialDataString == "LinearShear") initialData = initialDataLinearShear;
    if (InitialDataString == "Ring")        initialData = initialDataRing;
    // if (InitialDataString == "Zrake")       initialData = initialDataZrake;
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
    fs->setPlmTheta (double (user["plm"]));
    fs->setRiemannSolver (rs);
    fo->setConservationLaw (cl);
    mo->setMeshGeometry (mg);
    ss->setSourceTermsFunction (sourceTermsFunction);
    ss->setMeshOperator (mo);
    ss->setFieldOperator (fo);
    ss->setBoundaryCondition (bc);
    ss->setDisableFieldCT (true);
    ss->setIntercellFluxScheme (fs);
    ss->setViscousFluxFunction (makeViscousFluxCorrection (dx, dy));
    md->setVelocityIndex (cl->getIndexFor (ConservationLaw::VariableType::velocity));
    bc->setMeshGeometry (mg);
    bc->setConservationLaw (cl);

    if (LiveBinary)
    {
        ss->setRungeKuttaOrder (3);
        ss->setStarParticleDerivatives (starParticleDerivatives);
    }
    else
    {
        ss->setRungeKuttaOrder (2);
        ss->setStarParticleLocations (starParticleLocations);
    }

    auto L         = mo->linearCellDimension();
    auto V         = mo->measure (MeshLocation::cell);
    auto P         = md->getPrimitive();
    auto advance   = [&] (double dt) { return ss->advance (*md, status.simulationTime, dt); };
    auto condition = [&] () { return status.simulationTime < double (user["tfinal"]); };
    auto timestep  = [&] () { return timestepSize; };
    auto world     = MpiCommunicator::world();

    auto taskRecomputeDt = [&] (SimulationStatus, int rep)
    {
        const double dtGas  = fo->getCourantTimestep (P, L);
        const double dtStar = (LiveBinary > 0) ? HUGE_TIME_SCALE : dx / vStarMax;
        double localDt = double (user["cfl"]) * std::min (dtGas, dtStar);
        timestepSize = world.minimum (localDt);
    };

    auto taskCheckpoint = [&] (SimulationStatus, int rep)
    {
        writer->writeCheckpoint (rep, status, *cl, *md, *mg, *logger);

        if (int (user["tspurge"]))
        {
            tseries->clear();
        }
    };

    auto taskWriteGlobalStarParticleLocations = [&] (SimulationStatus, int rep)
    {
        for (std::size_t n = 0; n < std::min(md->starParticles.size(), std::size_t(8)); ++n)
        {
            globalStarParticleData[n] = md->starParticles[n];
        }
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

    scheduler->schedule (taskWriteGlobalStarParticleLocations, TaskScheduler::Recurrence (0.0, 0.0, 1), "write_global_star_particle");
    scheduler->schedule (taskCheckpoint, TaskScheduler::Recurrence (user["cpi"]), "checkpoint");
    scheduler->schedule (taskRecomputeDt, TaskScheduler::Recurrence (0.0, 0.0, 1), "compute_dt");
    scheduler->schedule (taskTimeSeries, TaskScheduler::Recurrence (user["tsi"]), "time_series");

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

        if (int (user["tspurge"]))
        {
            tseries->clear();
        }
    }
    else
    {
        md->starParticles = initialStarParticleLocations();
        taskWriteGlobalStarParticleLocations({}, 0);
        md->assignPrimitive (mo->generate (initialData, MeshLocation::cell));
    }
    md->applyBoundaryCondition (*bc);
    taskRecomputeDt (status, 0);

    return maraMainLoop (status, timestep, condition, advance, *scheduler, *logger);
}












































// ============================================================================
int BinaryTorqueStressCalculation::run (int argc, const char* argv[])
{
    for (int n = 1; n < argc; ++n)
    {
        processFile(argv[n]);
    }
    return 0;
}

void BinaryTorqueStressCalculation::processFile(std::string fname) const
{
    auto h5f  = H5::File (fname);
    auto user = h5f.getGroup ("user").readNamedValues();
    auto prim = h5f.getGroup ("primitive");
    auto stat = h5f.getGroup ("status");

    SofteningRadius   = user["SofteningRadius"];
    BetaBuffer        = user["BetaBuffer"];
    DomainRadius      = user["DomainRadius"];
    ViscousAlpha      = user["ViscousAlpha"];
    MachNumber        = user["MachNumber"];
    NumHoles          = user["NumHoles"];
    ScaleHeightFunc   = user["ScaleHeightFunc"];
    aBin              = user["aBin"];
    BinaryMassRatio   = user["BinaryMassRatio"];
    DiskMassRatio     = user["DiskMassRatio"];
    Eccentricity      = user["Eccentricity"];
    SinkRadius        = user["SinkRadius"];
    LocalViscousSink  = user["LocalViscousSink"];
    LdotEfficiency    = user["LdotEfficiency"];
    VacuumCleaner     = user["VacuumCleaner"];
    InitialDataString = std::string (user["InitialData"]);

    status.update (stat.readNamedValues());

    auto N  = int (user["N"]);
    auto cs = Shape {N, N, 1, 1, 1};
    auto bs = Shape {2, 2, 0, 0, 0};
    auto cl = std::make_shared<ThinDiskNewtonianHydro>();
    auto mg = std::make_shared<CartesianMeshGeometry>();
    auto rs = std::make_shared<HllcNewtonianHydroRiemannSolver>();
    auto fs = std::make_shared<MethodOfLinesPlm>();
    auto bc = std::make_shared<OutflowBoundaryCondition>();
    auto ss = std::make_shared<MethodOfLinesTVD>();
    auto mo = std::make_shared<MeshOperator>();
    auto fo = std::make_shared<FieldOperator>();
    auto md = std::make_shared<MeshData> (cs, bs, cl->getNumConserved());
    auto dx = mg->cellLength (0, 0, 0, 0);
    auto dy = mg->cellLength (0, 0, 0, 1);

    mg->setCellsShape (cs);
    mg->setLowerUpper ({-DomainRadius, -DomainRadius, 0.0}, {DomainRadius, DomainRadius, 1.0});
    bc->setMeshGeometry (mg);
    bc->setConservationLaw (cl);
    fs->setPlmTheta (user["plm"]);
    fs->setRiemannSolver (rs);
    fo->setConservationLaw (cl);
    mo->setMeshGeometry (mg);
    ss->setMeshOperator (mo);
    ss->setFieldOperator (fo);
    ss->setBoundaryCondition (bc);
    ss->setIntercellFluxScheme (fs);
    ss->setViscousFluxFunction (makeViscousFluxCorrection (dx, dy));

    auto P = prim.readArrays (cl->getPrimitiveNames(), 3);
    md->assignPrimitive (P);
    md->applyBoundaryCondition (*bc);

    fname.replace(fname.find(".h5"), 3, "-stress.h5");

    // Do something the following if this function is being used:
    assert(false); // it's disabled so it's not used incorrectly
    // for (std::size_t n = 0; n < std::min(md->starParticles.size(), std::size_t(8)); ++n)
    // {
    //     globalStarParticleData[n] = md->starParticles[n];
    // }

    auto h5res = H5::File (fname, "w");
    h5res.writeArray("viscous_flux", ss->computeViscousFluxes (*md, globalStarParticleData));
    h5res.writeArray("advective_flux", ss->computeAdvectiveFluxes (*md, globalStarParticleData));
}
