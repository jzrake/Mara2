#include <cmath>
#include "EulerEquation.hpp"

// Indexes to primitive quanitites P
#define RHO 0
#define V11 1
#define V22 2
#define V33 3
#define PRE 4

// Indexes to conserved quanitites U
#define DDD 0
#define S11 1
#define S22 2
#define S33 3
#define NRG 4




// ============================================================================
EulerEquation::EulerEquation() : gammaLawIndex (5./3)
{

}

ConservationLaw::State EulerEquation::fromConserved (const Request& request, const double* U) const
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

ConservationLaw::State EulerEquation::fromPrimitive (const Request& request, const double* P) const
{
    auto dAA = request.areaElement;
    const double gm0 = gammaLawIndex;
    const double gm1 = gammaLawIndex - 1.0;
    const double vv = P[V11] * P[V11] + P[V22] * P[V22] + P[V33] * P[V33];
    const double vn = P[V11] * dAA[0] + P[V22] * dAA[1] + P[V33] * dAA[2];
    const double cs = std::sqrt (gm0 * P[PRE] / P[RHO]);

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

int EulerEquation::getNumConserved() const
{
    return 5;
}
