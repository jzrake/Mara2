#include <cmath>
#include "ConservationLaws.hpp"

// Indexes to primitive quanitites P
#define RHO 0
#define V11 1
#define V22 2
#define V33 3
#define PRE 4
#define B11 5
#define B22 6
#define B33 7

// Indexes to conserved quanitites U
#define DDD 0
#define S11 1
#define S22 2
#define S33 3
#define NRG 4
#define H11 5
#define H22 6
#define H33 7




// ============================================================================
ScalarAdvection::ScalarAdvection (double waveSpeed) : waveSpeed (waveSpeed)
{

}

ConservationLaw::State ScalarAdvection::fromConserved (const Request& request, const double* U) const
{
    double u = U[0];
    double v = U[1];
    State S;
    S.P = {u, v};
    S.U = {u, v};
    S.A = {waveSpeed, waveSpeed};
    S.F = {waveSpeed * u, waveSpeed * v};
    S.L = Cow::Matrix (2, 2); // identity
    S.R = Cow::Matrix (2, 2);
    return S;
}

ConservationLaw::State ScalarAdvection::fromPrimitive (const Request& request, const double* P) const
{
    double u = P[0];
    double v = P[1];
    State S;
    S.P = {u, v};
    S.U = {u, v};
    S.A = {waveSpeed, waveSpeed};
    S.F = {waveSpeed * u, waveSpeed * v};
    S.L = Cow::Matrix (2, 2); // identity
    S.R = Cow::Matrix (2, 2);
    return S;
}

int ScalarAdvection::getNumConserved() const
{
    return 2;
}

std::string ScalarAdvection::getPrimitiveName (int fieldIndex) const
{
    switch (fieldIndex)
    {
        case 0: return "u";
        case 1: return "v";
        default: return "";
    }
}




// ============================================================================
NewtonianHydro::NewtonianHydro() : gammaLawIndex (5./3)
{

}

ConservationLaw::State NewtonianHydro::fromConserved (const Request& request, const double* U) const
{
    const double gm1 = gammaLawIndex - 1.0;
    const double pp = U[S11] * U[S11] + U[S22] * U[S22] + U[S33] * U[S33];
    double P[5];

    P[RHO] =  U[DDD];
    P[PRE] = (U[NRG] - 0.5 * pp / U[DDD]) * gm1;
    P[V11] =  U[S11] / U[DDD];
    P[V22] =  U[S22] / U[DDD];
    P[V33] =  U[S33] / U[DDD];

    return fromPrimitive (request, P);
}

ConservationLaw::State NewtonianHydro::fromPrimitive (const Request& request, const double* P) const
{
    const auto dAA = request.areaElement;
    const double gm0 = gammaLawIndex;
    const double gm1 = gammaLawIndex - 1.0;
    const double cs = std::sqrt (gm0 * P[PRE] / P[RHO]);
    const double vv = P[V11] * P[V11] + P[V22] * P[V22] + P[V33] * P[V33];
    const double vn = P[V11] * dAA[0] + P[V22] * dAA[1] + P[V33] * dAA[2];

    auto S = State();
    S.P.resize (5);
    S.U.resize (5);
    S.F.resize (5);
    S.A.resize (5);

    S.P[RHO] = P[RHO];
    S.P[V11] = P[V11];
    S.P[V22] = P[V22];
    S.P[V33] = P[V33];
    S.P[PRE] = P[PRE];

    S.U[DDD] = P[RHO];
    S.U[S11] = P[RHO] * P[V11];
    S.U[S22] = P[RHO] * P[V22];
    S.U[S33] = P[RHO] * P[V33];
    S.U[NRG] = P[RHO] * 0.5 * vv + P[PRE] / gm1;

    S.F[DDD] =  S.U[DDD] * vn;
    S.F[S11] =  S.U[S11] * vn + P[PRE] * dAA[0];
    S.F[S22] =  S.U[S22] * vn + P[PRE] * dAA[1];
    S.F[S33] =  S.U[S33] * vn + P[PRE] * dAA[2];
    S.F[NRG] = (S.U[NRG] + P[PRE]) * vn;

    S.A[0] = vn - cs;
    S.A[1] = vn;
    S.A[2] = vn;
    S.A[3] = vn;
    S.A[4] = vn + cs;

    return S;
}

int NewtonianHydro::getNumConserved() const
{
    return 5;
}

std::string NewtonianHydro::getPrimitiveName (int fieldIndex) const
{
    switch (fieldIndex)
    {
        case RHO: return "density";
        case S11: return "velocity1";
        case S22: return "velocity2";
        case S33: return "velocity3";
        case PRE: return "pressure";
        default: return "";
    }
}




// ============================================================================
NewtonianMHD::NewtonianMHD() : gammaLawIndex (5./3)
{

}

ConservationLaw::State NewtonianMHD::fromConserved (const Request& request, const double* U) const
{
    const double gm1 = gammaLawIndex - 1.0;
    const double pp = U[S11] * U[S11] + U[S22] * U[S22] + U[S33] * U[S33];
    const double BB = U[H11] * U[H11] + U[H22] * U[H22] + U[H33] * U[H33];
    double P[8];

    P[RHO] =  U[DDD];
    P[PRE] = (U[NRG] - 0.5 * pp / U[DDD] - 0.5 * BB) * gm1;
    P[V11] =  U[S11] / U[DDD];
    P[V22] =  U[S22] / U[DDD];
    P[V33] =  U[S33] / U[DDD];
    P[B11] =  U[H11];
    P[B22] =  U[H22];
    P[B33] =  U[H33];

    return fromPrimitive (request, P);
}

ConservationLaw::State NewtonianMHD::fromPrimitive (const Request& request, const double* P) const
{
    const auto dAA = request.areaElement;
    const double gm0 = gammaLawIndex;
    const double gm1 = gammaLawIndex - 1.0;
    const double cs = std::sqrt (gm0 * P[PRE] / P[RHO]);
    const double vv = P[V11] * P[V11] + P[V22] * P[V22] + P[V33] * P[V33];
    const double BB = P[B11] * P[B11] + P[B22] * P[B22] + P[B33] * P[B33];
    const double vn = P[V11] * dAA[0] + P[V22] * dAA[1] + P[V33] * dAA[2];
    const double Bn = P[B11] * dAA[0] + P[B22] * dAA[1] + P[B33] * dAA[2];

    auto S = State();
    S.P.resize (8);
    S.U.resize (8);
    S.F.resize (8);
    S.A.resize (8);

    S.P[RHO] = P[RHO];
    S.P[V11] = P[V11];
    S.P[V22] = P[V22];
    S.P[V33] = P[V33];
    S.P[PRE] = P[PRE];
    S.P[B11] = P[B11];
    S.P[B22] = P[B22];
    S.P[B33] = P[B33];

    S.U[DDD] = P[RHO];
    S.U[S11] = P[RHO] * P[V11];
    S.U[S22] = P[RHO] * P[V22];
    S.U[S33] = P[RHO] * P[V33];
    S.U[NRG] = P[RHO] * 0.5 * vv + 0.5 * BB + P[PRE] / gm1;
    S.U[H11] = P[B11];
    S.U[H22] = P[B22];
    S.U[H33] = P[B33];

    S.F[DDD] =  S.U[DDD] * vn;
    S.F[S11] =  S.U[S11] * vn + P[PRE] * dAA[0];
    S.F[S22] =  S.U[S22] * vn + P[PRE] * dAA[1];
    S.F[S33] =  S.U[S33] * vn + P[PRE] * dAA[2];
    S.F[NRG] = (S.U[NRG] + P[PRE]) * vn;
    S.F[H11] =  P[B11] * vn - Bn * P[V11];
    S.F[H22] =  P[B22] * vn - Bn * P[V22];
    S.F[H33] =  P[B33] * vn - Bn * P[V33];

    const double Bn2 = Bn * Bn;
    const double cs2 = cs * cs;
    const double ca2 = BB / P[RHO]; /* Alfven */
    const double cw4 = (cs2 + ca2) * (cs2 + ca2);
    const double cF2 = 0.5 * (cs2 + ca2 + std::sqrt (cw4 - 4 * cs2 * Bn2 / P[RHO])); /* fast */
    const double cS2 = 0.5 * (cs2 + ca2 - std::sqrt (cw4 - 4 * cs2 * Bn2 / P[RHO])); /* slow */

    S.A[0] = vn - std::sqrt (cF2);
    S.A[1] = vn - std::sqrt (ca2);
    S.A[2] = vn - std::sqrt (cS2);
    S.A[3] = vn;
    S.A[4] = vn;
    S.A[5] = vn + std::sqrt (cS2);
    S.A[6] = vn + std::sqrt (ca2);
    S.A[7] = vn + std::sqrt (cF2);

    return S;
}

int NewtonianMHD::getNumConserved() const
{
    return 8;
}

std::string NewtonianMHD::getPrimitiveName (int fieldIndex) const
{
    switch (fieldIndex)
    {
        case RHO: return "density";
        case S11: return "velocity1";
        case S22: return "velocity2";
        case S33: return "velocity3";
        case PRE: return "pressure";
        case B11: return "magnetic1";
        case B22: return "magnetic2";
        case B33: return "magnetic3";
        default: return "";
    }
}
