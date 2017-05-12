#include <cmath>
#include "Reconstruction.hpp"

#define MIN2(a, b) ((a) < (b) ? a : b)
#define MIN3(a, b, c) ((a) < (b) ? MIN2(a, c) : MIN2(b, c))
#define MAX2(a, b) ((a) > (b) ? a : b)
#define MAX3(a, b, c) ((a) > (b) ? MAX2(a, c) : MAX2(b, c))
#define MIN3ABS(a, b, c) MIN3(std::fabs(a), std::fabs(b), std::fabs(c));
#define SQU(x) ((x) * (x))
#define SGN(x) (((x) > 0) - ((x) < 0))




// ============================================================================
Reconstruction::Reconstruction()
{
    plmTheta = 2.0;
    shenzha10A = 50.0;
    modeIS = OriginalJiangShu96;
}

double Reconstruction::reconstruct (const double* v, enum Operation type) const
{
    switch (type)
    {
        case WENO5_FD_C2L: return weno5 (v, CeesC2L_FD, DeesC2L_FD);
        case WENO5_FD_C2R: return weno5 (v, CeesC2R_FD, DeesC2R_FD);
        case WENO5_FV_C2L: return weno5 (v, CeesC2L_FV, DeesC2L_FV);
        case WENO5_FV_C2R: return weno5 (v, CeesC2R_FV, DeesC2R_FV);
        case WENO5_FV_A2C: return weno5 (v, CeesA2C_FV, DeesA2C_FV);
        case WENO5_FV_C2A: return weno5 (v, CeesC2A_FV, DeesC2A_FV);
        case PLM_C2L: return plm (v, -1.0);
        case PLM_C2R: return plm (v, +1.0);
        default: return 0.0;
    }
}

void Reconstruction::setSmoothnessIndicator (enum SmoothnessIndicator IS)
{
    modeIS = IS;
}

void Reconstruction::setPlmTheta (double theta)
{
    plmTheta = theta;
}

void Reconstruction::setShenzha10A (double A)
{
    shenzha10A = A;
}

double Reconstruction::minmod (double ul, double u0, double ur) const
{
    const double a = plmTheta * (u0 - ul);
    const double b =      0.5 * (ur - ul);
    const double c = plmTheta * (ur - u0);
    return 0.25 * std::fabs (SGN(a) + SGN(b)) * (SGN(a) + SGN(c)) * MIN3ABS(a, b, c);
}

double Reconstruction::plm (const double* v, double sgn) const
{
    return v[0] + sgn * 0.5 * minmod (v[-1], v[0], v[1]);
}

double Reconstruction::weno5 (const double* v, const double c[3][3], const double d[3]) const
{
    double eps = 1e-6;
    double eps_prime = 1e-6;
    double w[3];

    const double vs[3] =
    {
        c[0][0] * v[+0] + c[0][1] * v[+1] + c[0][2] * v[2],
        c[1][0] * v[-1] + c[1][1] * v[+0] + c[1][2] * v[1],
        c[2][0] * v[-2] + c[2][1] * v[-1] + c[2][2] * v[0],
    };
    double B[3] = // smoothness indicators
    {
        (13./12.) * SQU(1 * v[+0] - 2 * v[+1] + 1 * v[+2]) +
        ( 1./ 4.) * SQU(3 * v[+0] - 4 * v[+1] + 1 * v[+2]),
        (13./12.) * SQU(1 * v[-1] - 2 * v[+0] + 1 * v[+1]) +
        ( 1./ 4.) * SQU(1 * v[-1] - 0 * v[+0] - 1 * v[+1]),
        (13./12.) * SQU(1 * v[-2] - 2 * v[-1] + 1 * v[+0]) +
        ( 1./ 4.) * SQU(1 * v[-2] - 4 * v[-1] + 3 * v[+0])
    };

    if (modeIS == ImprovedBorges08)
    {
        eps = eps_prime = 1e-14; // Borges uses 1e-40, but has Matlab
        const double tau5 = fabs(B[0] - B[2]);

        // Calculate weights with new smoothness indicators accoding to Borges
        w[0] = d[0] * (1.0 + (tau5 / (B[0] + eps)));
        w[1] = d[1] * (1.0 + (tau5 / (B[1] + eps)));
        w[2] = d[2] * (1.0 + (tau5 / (B[2] + eps)));
    }
    else if (modeIS == ImprovedShenZha10)
    {
        eps = 1e-6;
        eps_prime = 1e-10;
        const double A = shenzha10A; // [0 (less aggressive) -> ~100 (more aggressive)]
        const double minB = MIN3(B[0], B[1], B[2]);
        const double maxB = MAX3(B[0], B[1], B[2]);
        const double R0 = minB / (maxB + eps_prime);
        B[0] = R0 * A * minB + B[0];
        B[1] = R0 * A * minB + B[1];
        B[2] = R0 * A * minB + B[2];
        w[0] = d[0] / SQU(eps_prime + B[0]);
        w[1] = d[1] / SQU(eps_prime + B[1]);
        w[2] = d[2] / SQU(eps_prime + B[2]);
    }
    else // Use OriginalJiangShu96
    {
        eps = eps_prime = 1e-6; // recommended value by Jiang and Shu
        w[0] = d[0] / SQU(eps_prime + B[0]);
        w[1] = d[1] / SQU(eps_prime + B[1]);
        w[2] = d[2] / SQU(eps_prime + B[2]);
    }

    const double wtot = w[0] + w[1] + w[2];
    return (w[0] * vs[0] + w[1] * vs[1] + w[2] * vs[2]) / wtot;
}

const double Reconstruction::CeesA2C_FV[3][3] = {
    {23./24.,  1./12.,  -1./24.},
    {-1./24., 13./12.,  -1./24.},
    {-1./24.,  1./12.,  23./24.} };

const double Reconstruction::CeesC2A_FV[3][3] = {
    {25./24., -1./12.,   1./24.},
    { 1./24., 11./12.,   1./24.},
    { 1./24., -1./12.,  25./24.} };

const double Reconstruction::CeesC2L_FV[3][3] = {
    {15./8., -5./4.,  3./8.},
    { 3./8.,  3./4., -1./8.},
    {-1./8.,  3./4.,  3./8.} };

const double Reconstruction::CeesC2R_FV[3][3] = {
    { 3./8., 3./4.,  -1./8.},
    {-1./8., 3./4.,   3./8.},
    { 3./8.,-5./4.,  15./8.} };

const double Reconstruction::CeesC2L_FD[3][3] = {
    {11./6., -7./6.,  1./3. },
    { 1./3.,  5./6., -1./6. },
    {-1./6.,  5./6.,  1./3. } };

const double Reconstruction::CeesC2R_FD[3][3] = {
    { 1./3.,  5./6., -1./6. },
    {-1./6.,  5./6.,  1./3. },
    { 1./3., -7./6., 11./6. } };

const double Reconstruction::DeesA2C_FV[3] = {  -9./ 80.,  49./ 40.,  -9./ 80. };
const double Reconstruction::DeesC2A_FV[3] = { -17./240., 137./120., -17./240. };
const double Reconstruction::DeesC2L_FV[3] = {   1./ 16.,   5./  8.,   5./ 16. };
const double Reconstruction::DeesC2R_FV[3] = {   5./ 16.,   5./  8.,   1./ 16. };
const double Reconstruction::DeesC2L_FD[3] = { 0.1, 0.6, 0.3 };
const double Reconstruction::DeesC2R_FD[3] = { 0.3, 0.6, 0.1 };
