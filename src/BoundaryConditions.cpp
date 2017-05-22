#include <cmath>
#include "BoundaryConditions.hpp"




// ============================================================================
void PeriodicBoundaryCondition::apply (Cow::Array& P, int numGuard) const
{
    if (P.size(0) > 1) applyToAxis (P, numGuard, 0);
    if (P.size(1) > 1) applyToAxis (P, numGuard, 1);
    if (P.size(2) > 1) applyToAxis (P, numGuard, 2);
}

void PeriodicBoundaryCondition::applyToAxis (Cow::Array& P, int numGuard, int axis) const
{
    if (P.size (axis) == 1)
    {
        throw std::runtime_error("Attempt to apply boundary "
            "condition on flattened axis " + std::to_string (axis));
    }

    const int ng = numGuard;
    const int ni = P.size(0) - 2 * ng;
    const int nj = P.size(1) - 2 * ng;
    const int nk = P.size(2) - 2 * ng;
    const int nq = P.size(3);

    if (axis == 0)
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

    if (axis == 1)
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

    if (axis == 2)
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
void OutflowBoundaryCondition::apply (Cow::Array& P, int numGuard) const
{
    if (P.size(0) > 1) applyToAxis (P, numGuard, 0);
    if (P.size(1) > 1) applyToAxis (P, numGuard, 1);
    if (P.size(2) > 1) applyToAxis (P, numGuard, 2);
}

void OutflowBoundaryCondition::applyToAxis (Cow::Array& P, int numGuard, int axis) const
{
    if (P.size (axis) == 1)
    {
        throw std::runtime_error("Attempt to apply boundary "
            "condition on flattened axis " + std::to_string (axis));
    }

    const int ng = numGuard;
    const int ni = P.size(0) - 2 * ng;
    const int nj = P.size(1) - 2 * ng;
    const int nk = P.size(2) - 2 * ng;
    const int nq = P.size(3);

    if (axis == 0)
    {
        for (int j = 0; j < P.size(1); ++j)
        {
            for (int k = 0; k < P.size(2); ++k)
            {
                for (int i = 0; i < ng; ++i)
                {
                    for (int q = 0; q < nq; ++q)
                    {
                        P (i,           j, k, q) = P (ng,          j, k, q);
                        P (i + ni + ng, j, k, q) = P (ng + ni - 1, j, k, q);
                    }
                }
            }
        }
    }

    if (axis == 1)
    {
        for (int k = 0; k < P.size(2); ++k)
        {
            for (int i = 0; i < P.size(0); ++i)
            {
                for (int j = 0; j < ng; ++j)
                {
                    for (int q = 0; q < nq; ++q)
                    {
                        P (i, j,           k, q) = P (i, ng,          k, q);
                        P (i, j + nj + ng, k, q) = P (i, nj + ng - 1, k, q);
                    }
                }
            }
        }
    }

    if (axis == 2)
    {
        for (int i = 0; i < P.size(0); ++i)
        {
            for (int j = 0; j < P.size(1); ++j)
            {
                for (int k = 0; k < ng; ++k)
                {
                    for (int q = 0; q < nq; ++q)
                    {
                        P (i, j, k,           q) = P (i, j, ng,          q);
                        P (i, j, k + nk + ng, q) = P (i, j, ng + nk - 1, q);
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
    auto periodic = PeriodicBoundaryCondition();
    auto outflow = OutflowBoundaryCondition();

    if (P.size(0) > 1) periodic.applyToAxis (P, numGuard, 0);
    if (P.size(1) > 1) periodic.applyToAxis (P, numGuard, 1);
    if (P.size(2) > 1) outflow.applyToAxis (P, numGuard, 2);

    const int V11 = 1;
    const int V22 = 2;
    const int ng = numGuard;
    const int ni = P.size(0) - 2 * ng;
    const int nj = P.size(1) - 2 * ng;
    const int nk = P.size(2) - 2 * ng;

    for (int i = 0; i < P.size(0); ++i)
    {
        for (int j = 0; j < P.size(1); ++j)
        {
            const double x = (i - ng + 0.5) / ni - 0.5;
            const double y = (j - ng + 0.5) / nj - 0.5;
            const double vx = 0.1 * std::sin (2 * M_PI * y);
            const double vy = 0.1 * std::cos (2 * M_PI * x);

            for (int k = 0; k < ng; ++k)
            {
                P (i, j, k,           V11) =  vx;
                P (i, j, k + nk + ng, V11) = -vx;
                P (i, j, k,           V22) =  vy;
                P (i, j, k + nk + ng, V22) = -vy;
            }
        }
    }
}
