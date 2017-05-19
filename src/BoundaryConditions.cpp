#include <cmath>
#include "BoundaryConditions.hpp"




// ============================================================================
void PeriodicBoundaryCondition::apply (Cow::Array& P, int numGuard) const
{
    const int ng = numGuard;
    const int ni = P.size(0) - 2 * ng;
    const int nj = P.size(1) - 2 * ng;
    const int nk = P.size(2) - 2 * ng;
    const int nq = P.size(3);

    if (P.size(0) > 1)
    {
        for (int j = 0; j < P.size(1); ++j)
        {
            for (int k = 0; k < P.size(2); ++k)
            {
                for (int i = 0; i < ng; ++i)
                {
                    for (int q = 0; q < nq; ++q)
                    {
                        P (i,           j, k, q) = P (i + ni, j, k, q);
                        P (i + ng + ni, j, k, q) = P (i + ng, j, k, q);
                    }
                }
            }
        }
    }

    if (P.size(1) > 1)
    {
        for (int k = 0; k < P.size(2); ++k)
        {
            for (int i = 0; i < P.size(0); ++i)
            {
                for (int j = 0; j < ng; ++j)
                {
                    for (int q = 0; q < nq; ++q)
                    {
                        P (i, j,           k, q) = P (i, j + nj, k, q);
                        P (i, j + ng + nj, k, q) = P (i, j + ng, k, q);
                    }
                }
            }
        }
    }

    if (P.size(2) > 1)
    {
        for (int i = 0; i < P.size(0); ++i)
        {
            for (int j = 0; j < P.size(1); ++j)
            {
                for (int k = 0; k < ng; ++k)
                {
                    for (int q = 0; q < nq; ++q)
                    {
                        P (i, j, k,           q) = P (i, j, k + nk, q);
                        P (i, j, k + ng + nk, q) = P (i, j, k + ng, q);
                    }
                }
            }
        }
    }
}




// ============================================================================
void PlanarPipeFlow::apply (Cow::Array& P, int numGuard) const
{
    const int RHO = 0;
    const int V11 = 1;
    const int V22 = 2;
    const int V33 = 3;
    const int PRE = 4;
    const int ng = numGuard;
    const int ni = P.size(0) - 2 * ng;
    const int nj = P.size(1) - 2 * ng;
    const double inflowVelocity = 0.5;

    for (int j = 0; j < P.size(1); ++j)
    {
        for (int i = 0; i < ng; ++i)
        {
            // Inflow on the left side of the domain
            P (i, j, 0, RHO) = P (ng, j, 0, RHO);
            P (i, j, 0, V11) = inflowVelocity;
            P (i, j, 0, V22) = P (ng, j, 0, V22);
            P (i, j, 0, V33) = P (ng, j, 0, V33);
            P (i, j, 0, PRE) = P (ng, j, 0, PRE);

            // Outflow on the right side of the domain
            P (ni + ng + i, j, 0, RHO) = P (ni + ng - 1, j, 0, RHO);
            P (ni + ng + i, j, 0, V11) = P (ni + ng - 1, j, 0, V11);
            P (ni + ng + i, j, 0, V22) = P (ni + ng - 1, j, 0, V22);
            P (ni + ng + i, j, 0, V33) = P (ni + ng - 1, j, 0, V33);
            P (ni + ng + i, j, 0, PRE) = P (ni + ng - 1, j, 0, PRE);
        }
    }

    // Reflecting BC's on the top and bottom walls
    if (P.size (1) > 1)
    {
        for (int i = 0; i < P.size(0); ++i)
        {
            for (int j = 0; j < ng; ++j)
            {
                // Bottom wall
                P (i, ng - (j + 1), 0, RHO) = P (i, ng + j, 0, RHO);
                P (i, ng - (j + 1), 0, V11) = 0.0;
                P (i, ng - (j + 1), 0, V22) = P (i, ng + j, 0, V22) * -1.0;
                P (i, ng - (j + 1), 0, V33) = 0.0;
                P (i, ng - (j + 1), 0, PRE) = P (i, ng + j, 0, PRE);

                // Top wall
                P (i, nj + ng - 1 + (j + 1), 0, RHO) = P (i, nj + ng - 1 - j, 0, RHO);
                P (i, nj + ng - 1 + (j + 1), 0, V11) = 0.0;
                P (i, nj + ng - 1 + (j + 1), 0, V22) = P (i, nj + ng - 1 - j, 0, V22) * -1.0;
                P (i, nj + ng - 1 + (j + 1), 0, V33) = 0.0;
                P (i, nj + ng - 1 + (j + 1), 0, PRE) = P (i, nj + ng - 1 - j, 0, PRE);
            }
        }
    }
}




// ============================================================================
void DrivenMHDBoundary::apply (Cow::Array& P, int numGuard) const
{
    PeriodicBoundaryCondition periodic;
    periodic.apply (P, numGuard);

    const int RHO = 0;
    const int V11 = 1;
    const int V22 = 2;
    const int V33 = 3;
    const int PRE = 4;
    const int B11 = 5;
    const int B22 = 6;
    const int B33 = 7;
    const int ng = numGuard;
    const int ni = P.size(0) - 2 * ng;
    const int nj = P.size(1) - 2 * ng;
    const int nk = P.size(2) - 2 * ng;
    const double v0 = 0.1;

    for (int i = ng; i < ni + ng; ++i)
    {
        for (int j = ng; j < nj + ng; ++j)
        {
            const double x = (i + 0.5 - ng - 0.5 * ni) / ni;
            const double y = (j + 0.5 - ng - 0.5 * nj) / nj;

            for (int k = 0; k < ng; ++k)
            {
                P (i, j, ng - (k + 1), RHO) = P (i, j, ng + k, RHO);
                P (i, j, ng - (k + 1), V11) = v0 * std::sin (y);
                P (i, j, ng - (k + 1), V22) = v0 * std::cos (x);
                P (i, j, ng - (k + 1), V33) = P (i, j, ng + k, V33) * -1.0;
                P (i, j, ng - (k + 1), PRE) = P (i, j, ng + k, PRE);
                P (i, j, ng - (k + 1), B11) = P (i, j, ng, B11);
                P (i, j, ng - (k + 1), B22) = P (i, j, ng, B22);
                P (i, j, ng - (k + 1), B33) = P (i, j, ng, B33);

                P (i, j, nk + ng - 1 + (k + 1), RHO) = P (i, j, nk + ng - 1 - k, RHO);
                P (i, j, nk + ng - 1 + (k + 1), V11) = v0 * std::sin (y) * -1.0;
                P (i, j, nk + ng - 1 + (k + 1), V22) = v0 * std::cos (x) * -1.0;
                P (i, j, nk + ng - 1 + (k + 1), V33) = P (i, j, nk + ng - 1 - k, V33) * -1.0;
                P (i, j, nk + ng - 1 + (k + 1), PRE) = P (i, j, nk + ng - 1 - k, PRE);
                P (i, j, ng + (k + 1), B11) = P (i, j, ng, B11);
                P (i, j, ng + (k + 1), B22) = P (i, j, ng, B22);
                P (i, j, ng + (k + 1), B33) = P (i, j, ng, B33);
            }
        }
    }
}
