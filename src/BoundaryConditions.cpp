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
// #include <cmath>

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
                P (i, ng - (j + 1), 0, RHO) = P (i, ng + (j + 1), 0, RHO);
                P (i, ng - (j + 1), 0, V11) = 0.0;
                P (i, ng - (j + 1), 0, V22) = P (i, ng + (j + 1), 0, V22) * -1.0;
                P (i, ng - (j + 1), 0, V33) = 0.0;
                P (i, ng - (j + 1), 0, PRE) = P (i, ng + (j + 1), 0, PRE);

                P (i, nj + ng - 1 + (j + 1), 0, RHO) = P (i, nj + ng - 1 - (j + 1), 0, RHO);
                P (i, nj + ng - 1 + (j + 1), 0, V11) = 0.0;
                P (i, nj + ng - 1 + (j + 1), 0, V22) = P (i, nj + ng - 1 - (j + 1), 0, V22) * -1.0;
                P (i, nj + ng - 1 + (j + 1), 0, V33) = 0.0;
                P (i, nj + ng - 1 + (j + 1), 0, PRE) = P (i, nj + ng - 1 - (j + 1), 0, PRE);
            }
        }
    }
}




// ============================================================================
void DrivenMHDBoundary::apply (Cow::Array& P, int numGuard) const
{
    PeriodicBoundaryCondition periodic;
    periodic.apply (P, numGuard);
}
