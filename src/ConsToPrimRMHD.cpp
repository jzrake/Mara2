/*------------------------------------------------------------------------------
 * FILE: ConsToPrimRMHD.cpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 * REFERENCES:
 *
 * Anton, L., Zanotti, O., Miralles, J. A., Martı, J. M., Ibanez, J. M., Font,
 * J. A., & Pons, J. A. 2006, The Astrophysical Journal, 637, 296
 *
 * Noble, S. C., Gammie, C. F., McKinney, J. C., & Zanna, L. D.  2006, The
 * Astrophysical Journal, 641, 626
 *
 *
 * DESCRIPTION:
 *
 * This module solves for the primitive state (rho, pre, v, B) given the
 * conserved state (ddd, tau, S, B). In other words it solves for the 5 unknowns
 * rho, pre, vx, vy, vz. An adiabatic equation of state is assumed throughout,
 * with the index being provided by srmhd_c2p_set_gamma(). The solve functions
 * promise not to modify the pointer to result primitives unless the execution
 * is successful.
 *
 * ------------------------------------------------------------------------------
 */


#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cmath>
#include "ConsToPrimRMHD.hpp"



enum { ddd, Sx, Sy, Sz, tau, Bx, By, Bz }; // Conserved
enum { rho, vx, vy, vz, pre };             // Primitive



ConsToPrimRMHD::ConsToPrimRMHD()
{
    MaxIterations = 250;
    Tolerance = 1e-12;
    bigZ = 1e20;
    bigW = 1e12;
    smlZ = 0.0;
    smlW = 1.0;
    Iterations = 0;
    AppliedPressureFloor = 0;
    AdiabaticGamma = 4. / 3;
    PressureFloor = -1.0;
}

int ConsToPrimRMHD::did_use_pressure_floor()
{
    return AppliedPressureFloor;
}

void ConsToPrimRMHD::set_gamma (double adiabatic_gamma)
{
    AdiabaticGamma = adiabatic_gamma;
    gamf = (AdiabaticGamma - 1.0) / AdiabaticGamma;
}

int ConsToPrimRMHD::get_iterations()
{
    return Iterations;
}

int ConsToPrimRMHD::set_pressure_floor (double pf)
{
    PressureFloor = pf;
    return pf < 0.0;
}

void ConsToPrimRMHD::new_state (const double *U)
{
    D    = U[ddd];
    Tau  = U[tau];
    S2   = U[Sx]*U[Sx] + U[Sy]*U[Sy] + U[Sz]*U[Sz];
    B2   = U[Bx]*U[Bx] + U[By]*U[By] + U[Bz]*U[Bz];
    BS   = U[Bx]*U[Sx] + U[By]*U[Sy] + U[Bz]*U[Sz];
    BS2  = BS * BS;
    std::memcpy (Cons, U, 8*sizeof(double));
}

void ConsToPrimRMHD::estimate_from_cons()
{
    Z_start = std::sqrt (S2 + D*D);
    W_start = Z_start / D;
}

void ConsToPrimRMHD::get_starting_prim (double *P)
{
    reconstruct_prim (Z_start, W_start, P);
}

void ConsToPrimRMHD::set_starting_prim (const double *P)
{
    const double V2 = P[vx]*P[vx] + P[vy]*P[vy] + P[vz]*P[vz];
    const double W2 = 1.0 / (1.0 - V2);
    const double e = P[pre] / (P[rho] * (AdiabaticGamma - 1.0));
    const double h = 1.0 + e + P[pre] / P[rho];

    Z_start = P[rho] * h * W2;
    W_start = std::sqrt (W2);
}

int ConsToPrimRMHD::reconstruct_prim (double Z, double W, double *Pout)
{
    double P[8]; // Place result into temporary prim state for now.
    const double b0 = BS * W / Z;

    P[rho] =  D/W;
    P[pre] = (D/W) * (Z/(D*W) - 1.0) * gamf;
    P[vx ] = (Cons[Sx] + b0 * Cons[Bx] / W) / (Z + B2);
    P[vy ] = (Cons[Sy] + b0 * Cons[By] / W) / (Z + B2);
    P[vz ] = (Cons[Sz] + b0 * Cons[Bz] / W) / (Z + B2);
    P[Bx ] =  Cons[Bx];
    P[By ] =  Cons[By];
    P[Bz ] =  Cons[Bz];

    if (P[pre] <= 0.0 && PressureFloor > 0.0)
    {
        P[pre] = P[rho] * PressureFloor;
        AppliedPressureFloor = 1;
    }
    else
    {
        AppliedPressureFloor = 0;
    }
    int error = check_prim(P);

    if (error)
    {
        return error;
    }

    // Don't actually modify the user's prim state until we're sure the state
    // is good.
    // ------------------------------------------------------------------------
    std::memcpy (Pout, P, 8 * sizeof (double));
    return SRMHD_C2P_SUCCESS;
}

int ConsToPrimRMHD::check_cons (const double *U)
{
    int i;

    if (U[ddd] < 0.0) return SRMHD_C2P_CONS_NEGATIVE_DENSITY;
    if (U[tau] < 0.0) return SRMHD_C2P_CONS_NEGATIVE_ENERGY;

    for (i=0; i<8; ++i)
    {
        if (! std::isfinite (U[i]))
        {
            return SRMHD_C2P_CONS_CONTAINS_NAN;
        }
    }
    return SRMHD_C2P_SUCCESS;
}

int ConsToPrimRMHD::check_prim (const double *P)
{
    int i;
    const double v2 = P[vx]*P[vx] + P[vy]*P[vy] + P[vz]*P[vz];

    if (v2   >=  1.0) return SRMHD_C2P_PRIM_SUPERLUMINAL;
    if (P[pre] < 0.0) return SRMHD_C2P_PRIM_NEGATIVE_PRESSURE;
    if (P[rho] < 0.0) return SRMHD_C2P_PRIM_NEGATIVE_RESTMASS;

    for (i=0; i<8; ++i)
    {
        if (! std::isfinite (P[i]))
        {
            return SRMHD_C2P_PRIM_CONTAINS_NAN;
        }
  }
  return SRMHD_C2P_SUCCESS;
}

int ConsToPrimRMHD::solve_anton2dzw (double *P)
{
    if (int bad_input = check_cons (Cons))
    {
        return bad_input;
    }

    // Starting values
    // ---------------------------------------------------------------------------
    Iterations = 0;

    double error = 1.0;
    double W = W_start;
    double Z = Z_start;

    while (error > Tolerance)
    {
        const double Z2 = Z*Z;
        const double Z3 = Z*Z2;
        const double W2 = W*W;
        const double W3 = W*W2;
        const double Pre = (D/W) * (Z/(D*W) - 1.0) * gamf;

        const double df0dZ = 2*(B2+Z)*(BS2*W2 + (W2-1)*Z3) / (W2*Z3);
        const double df0dW = 2*(B2+Z)*(B2+Z) / W3;
        const double df1dZ = 1.0 + BS2/Z3 - gamf/W2;
        const double df1dW = B2/W3 + (2*Z - D*W)/W3 * gamf;

        double f[2];
        double J[4], G[4];

        // Evaluation of the function, and its Jacobian
        // -------------------------------------------------------------------------
        f[0] = -S2  + (Z+B2)*(Z+B2)*(W2-1)/W2 - (2*Z+B2)*BS2/Z2;      // eqn (84)
        f[1] = -Tau +  Z+B2 - Pre - 0.5*B2/W2 -      0.5*BS2/Z2 - D;  // eqn (85)

        J[0] = df0dZ; J[1] = df0dW;
        J[2] = df1dZ; J[3] = df1dW;

        // G in the inverse Jacobian
        // -------------------------------------------------------------------------
        const double det = J[0]*J[3] - J[1]*J[2];
        G[0] =  J[3]/det; G[1] = -J[1]/det;
        G[2] = -J[2]/det; G[3] =  J[0]/det;      // G = J^{-1}

        const double dZ = -(G[0]*f[0] + G[1]*f[1]); // Matrix multiply, dx = -G . f
        const double dW = -(G[2]*f[0] + G[3]*f[1]);

        // Bracketing the root
        // -------------------------------------------------------------------------
        double Z_new = Z + dZ;
        double W_new = W + dW;

        Z_new = (Z_new > smlZ) ? Z_new : -Z_new;
        Z_new = (Z_new < bigZ) ? Z_new :  Z;

        W_new = (W_new > smlW) ? W_new : smlW;
        W_new = (W_new < bigW) ? W_new : bigW;

        Z = Z_new;
        W = W_new;

        // -------------------------------------------------------------------------
        error = std::fabs (dZ/Z) + std::fabs (dW/W);
        ++Iterations;

        if (Iterations == MaxIterations)
        {
            return SRMHD_C2P_MAXITER;
        }
    }
    return reconstruct_prim (Z, W, P);
}

int ConsToPrimRMHD::solve_noble1dw (double *P)
{
    if (int bad_input = check_cons (Cons))
    {
        return bad_input;
    }

    // Starting values
    // ---------------------------------------------------------------------------
    Iterations = 0;

    double error = 1.0;
    double Z = Z_start;
    double f, g;

    while (error > Tolerance)
    {
        const double Z2  = Z*Z;
        const double Z3  = Z*Z2;
        const double a   = S2*Z2 + BS2*(B2 + 2*Z);
        const double b   = (B2 + Z)*(B2 + Z)*Z2;
        const double ap  = 2*(S2*Z + BS2);          // da/dZ
        const double bp  = 2*Z*(B2 + Z)*(B2 + 2*Z); // db/dZ
        const double V2  = a / b;
        const double W2  = 1.0 / (1.0 - V2);
        const double W   = std::sqrt (W2);
        const double W3  = W*W2;
        const double Pre = (D/W) * (Z/(D*W) - 1.0) * gamf;

        const double dv2dZ    = (ap*b - bp*a) / (b*b); // (a'b - b'a) / b^2
        const double delPdelZ = gamf/W2;
        const double delPdelW = gamf * (D/W2 - 2*Z/W3);
        const double dWdv2    = 0.5*W3;
        const double dPdZ     = delPdelW * dWdv2 * dv2dZ + delPdelZ;

        f = Tau + D - 0.5*B2*(1+V2) + 0.5*BS2/Z2 - Z + Pre; // equation (29)
        g = -0.5*B2*dv2dZ - BS2/Z3 - 1.0 + dPdZ;

        const double dZ = -f/g;

        // -------------------------------------------------------------------------
        double Z_new = Z + dZ;

        Z_new = (Z_new > smlZ) ? Z_new : -Z_new;
        Z_new = (Z_new < bigZ) ? Z_new :  Z;

        Z = Z_new;

        error = std::fabs (dZ / Z);
        ++Iterations;

        if (Iterations == MaxIterations)
        {
            return SRMHD_C2P_MAXITER;
        }
    }

    // Recover the W value from the converged Z value.
    // -------------------------------------------------------------------------
    const double Z2  = Z*Z;
    const double a   = S2*Z2 + BS2*(B2 + 2*Z);
    const double b   = (B2 + Z)*(B2 + Z)*Z2;
    const double V2  = a / b;
    const double W2  = 1.0 / (1.0 - V2);
    const double W   = std::sqrt (W2);
    return reconstruct_prim (Z, W, P);
}

const char *ConsToPrimRMHD::get_error (int error)
{
  switch (error) {
      case SRMHD_C2P_SUCCESS:
      return "successful inversion from conserved to primitive";
      case SRMHD_C2P_CONS_CONTAINS_NAN:
      return "input conserved state contained nan's";
      case SRMHD_C2P_CONS_NEGATIVE_DENSITY:
      return "input conserved state has negative density";
      case SRMHD_C2P_CONS_NEGATIVE_ENERGY:
      return "input conserved state has negative energy";
      case SRMHD_C2P_PRIM_CONTAINS_NAN:
      return "derived primitive state contains nan's";
      case SRMHD_C2P_PRIM_NEGATIVE_PRESSURE:
      return "derived primitive state has negative pressure";
      case SRMHD_C2P_PRIM_NEGATIVE_RESTMASS:
      return "derived primitive state contains negative density";
      case SRMHD_C2P_PRIM_SUPERLUMINAL:
      return "derived primitive state contains superluminal velocity";
      case SRMHD_C2P_MAXITER:
      return "rootfinder hit maximum iterations";
      default:
      return "unkown error";
  }
}
