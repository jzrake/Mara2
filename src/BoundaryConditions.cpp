#include <cmath>
#include "BoundaryConditions.hpp"




// ============================================================================
void BoundaryCondition::applyToCellCenteredB (Cow::Array& B, int numGuard) const
{
    throw std::runtime_error ("BoundaryCondition::applyToCellCenteredB required but not implemented");
}

void BoundaryCondition::applyToGodunovFluxes (Cow::Array& F, int numGuard, int axis) const
{
    throw std::runtime_error ("BoundaryCondition::applyToGodunovFluxes required but not implemented");
}




// ============================================================================
void PeriodicBoundaryCondition::apply (Cow::Array& P, const ConservationLaw& law, int numGuard) const
{
    if (P.size(0) > 1) applyToAxis (P, numGuard, 0);
    if (P.size(1) > 1) applyToAxis (P, numGuard, 1);
    if (P.size(2) > 1) applyToAxis (P, numGuard, 2);
}

void PeriodicBoundaryCondition::applyToCellCenteredB (Cow::Array& B, int numGuard) const
{
    if (B.size(0) > 1) applyToAxis (B, numGuard, 0);
    if (B.size(1) > 1) applyToAxis (B, numGuard, 1);
    if (B.size(2) > 1) applyToAxis (B, numGuard, 2);
}

void PeriodicBoundaryCondition::applyToGodunovFluxes (Cow::Array& F, int numGuard, int axis) const
{
    applyToAxis (F, numGuard, axis);
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
void OutflowBoundaryCondition::apply (Cow::Array& P, const ConservationLaw& law, int numGuard) const
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
void ReflectingBoundaryCondition::apply (Cow::Array& P, const ConservationLaw& law, int numGuard) const
{
    if (P.size(0) > 1) applyToAxis (P, law, numGuard, 0);
    if (P.size(1) > 1) applyToAxis (P, law, numGuard, 1);
    if (P.size(2) > 1) applyToAxis (P, law, numGuard, 2);
}

void ReflectingBoundaryCondition::applyToAxis (Cow::Array& P, const ConservationLaw& law, int numGuard, int axis) const
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

    const int irho = law.getIndexFor (ConservationLaw::VariableType::density);
    const int ipre = law.getIndexFor (ConservationLaw::VariableType::pressure);
    const int ivel = law.getIndexFor (ConservationLaw::VariableType::velocity);
    const int imag = law.getIndexFor (ConservationLaw::VariableType::magnetic);

    if (axis == 0)
    {
        for (int j = 0; j < P.size(1); ++j)
        {
            for (int k = 0; k < P.size(2); ++k)
            {
                for (int i = 0; i < ng; ++i)
                {
                    const int iExtL = +i;
                    const int iIntL = -i + 2 * ng - 1;
                    const int iIntR = +i + ni;
                    const int iExtR = -i + 2 * ng - 1 + ni;

                    if (irho != -1) P (iExtL, j, k, irho + 0) = P (iIntL, j, k, irho);
                    if (ipre != -1) P (iExtL, j, k, ipre + 0) = P (iIntL, j, k, ipre);
                    if (ivel != -1) P (iExtL, j, k, ivel + 0) = P (iIntL, j, k, ivel + 0) * -1.0;
                    if (ivel != -1) P (iExtL, j, k, ivel + 1) = P (iIntL, j, k, ivel + 1);
                    if (ivel != -1) P (iExtL, j, k, ivel + 2) = P (iIntL, j, k, ivel + 2);
                    if (imag != -1) P (iExtL, j, k, imag + 0) = P (iIntL, j, k, imag + 0);
                    if (imag != -1) P (iExtL, j, k, imag + 1) = P (iIntL, j, k, imag + 1) * -1.0;
                    if (imag != -1) P (iExtL, j, k, imag + 2) = P (iIntL, j, k, imag + 2) * -1.0;

                    if (irho != -1) P (iExtR, j, k, irho + 0) = P (iIntR, j, k, irho);
                    if (ipre != -1) P (iExtR, j, k, ipre + 0) = P (iIntR, j, k, ipre);
                    if (ivel != -1) P (iExtR, j, k, ivel + 0) = P (iIntR, j, k, ivel + 0) * -1.0;
                    if (ivel != -1) P (iExtR, j, k, ivel + 1) = P (iIntR, j, k, ivel + 1);
                    if (ivel != -1) P (iExtR, j, k, ivel + 2) = P (iIntR, j, k, ivel + 2);
                    if (imag != -1) P (iExtR, j, k, imag + 0) = P (iIntR, j, k, imag + 0);
                    if (imag != -1) P (iExtR, j, k, imag + 1) = P (iIntR, j, k, imag + 1);
                    if (imag != -1) P (iExtR, j, k, imag + 2) = P (iIntR, j, k, imag + 2);
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
                    const int jExtL = +j;
                    const int jIntL = -j + 2 * ng - 1;
                    const int jIntR = +j + nj;
                    const int jExtR = -j + 2 * ng - 1 + nj;

                    if (irho != -1) P (i, jExtL, k, irho + 0) = P (i, jIntL, k, irho);
                    if (ipre != -1) P (i, jExtL, k, ipre + 0) = P (i, jIntL, k, ipre);
                    if (ivel != -1) P (i, jExtL, k, ivel + 0) = P (i, jIntL, k, ivel + 0);
                    if (ivel != -1) P (i, jExtL, k, ivel + 1) = P (i, jIntL, k, ivel + 1) * -1.0;
                    if (ivel != -1) P (i, jExtL, k, ivel + 2) = P (i, jIntL, k, ivel + 2);
                    if (imag != -1) P (i, jExtL, k, imag + 0) = P (i, jIntL, k, imag + 0);
                    if (imag != -1) P (i, jExtL, k, imag + 1) = P (i, jIntL, k, imag + 1);
                    if (imag != -1) P (i, jExtL, k, imag + 2) = P (i, jIntL, k, imag + 2);

                    if (irho != -1) P (i, jExtR, k, irho + 0) = P (i, jIntR, k, irho);
                    if (ipre != -1) P (i, jExtR, k, ipre + 0) = P (i, jIntR, k, ipre);
                    if (ivel != -1) P (i, jExtR, k, ivel + 0) = P (i, jIntR, k, ivel + 0);
                    if (ivel != -1) P (i, jExtR, k, ivel + 1) = P (i, jIntR, k, ivel + 1) * -1.0;
                    if (ivel != -1) P (i, jExtR, k, ivel + 2) = P (i, jIntR, k, ivel + 2);
                    if (imag != -1) P (i, jExtR, k, imag + 0) = P (i, jIntR, k, imag + 0);
                    if (imag != -1) P (i, jExtR, k, imag + 1) = P (i, jIntR, k, imag + 1);
                    if (imag != -1) P (i, jExtR, k, imag + 2) = P (i, jIntR, k, imag + 2);
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
                    const int kExtL = +k;
                    const int kIntL = -k + 2 * ng - 1;
                    const int kIntR = +k + nk;
                    const int kExtR = -k + 2 * ng - 1 + nk;

                    if (irho != -1) P (i, j, kExtL, irho + 0) = P (i, j, kIntL, irho);
                    if (ipre != -1) P (i, j, kExtL, ipre + 0) = P (i, j, kIntL, ipre);
                    if (ivel != -1) P (i, j, kExtL, ivel + 0) = P (i, j, kIntL, ivel + 0);
                    if (ivel != -1) P (i, j, kExtL, ivel + 1) = P (i, j, kIntL, ivel + 1);
                    if (ivel != -1) P (i, j, kExtL, ivel + 2) = P (i, j, kIntL, ivel + 2) * -1.0;
                    if (imag != -1) P (i, j, kExtL, imag + 0) = P (i, j, kIntL, imag + 0);
                    if (imag != -1) P (i, j, kExtL, imag + 1) = P (i, j, kIntL, imag + 1);
                    if (imag != -1) P (i, j, kExtL, imag + 2) = P (i, j, kIntL, imag + 2);

                    if (irho != -1) P (i, j, kExtR, irho + 0) = P (i, j, kIntR, irho);
                    if (ipre != -1) P (i, j, kExtR, ipre + 0) = P (i, j, kIntR, ipre);
                    if (ivel != -1) P (i, j, kExtR, ivel + 0) = P (i, j, kIntR, ivel + 0);
                    if (ivel != -1) P (i, j, kExtR, ivel + 1) = P (i, j, kIntR, ivel + 1);
                    if (ivel != -1) P (i, j, kExtR, ivel + 2) = P (i, j, kIntR, ivel + 2) * -1.0;
                    if (imag != -1) P (i, j, kExtR, imag + 0) = P (i, j, kIntR, imag + 0);
                    if (imag != -1) P (i, j, kExtR, imag + 1) = P (i, j, kIntR, imag + 1);
                    if (imag != -1) P (i, j, kExtR, imag + 2) = P (i, j, kIntR, imag + 2);
                }
            }
        }
    }
}




// ============================================================================
void PlanarPipeFlow::apply (Cow::Array& P, const ConservationLaw& law, int numGuard) const
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
void DrivenMHDBoundary::apply (Cow::Array& P, const ConservationLaw& law, int numGuard) const
{
    auto periodic = PeriodicBoundaryCondition();
    auto vertical = ReflectingBoundaryCondition();

    if (P.size(0) > 1) periodic.applyToAxis (P, numGuard, 0);
    if (P.size(1) > 1) periodic.applyToAxis (P, numGuard, 1);
    if (P.size(2) > 1) vertical.applyToAxis (P, law, numGuard, 2);

    const int ivel = law.getIndexFor (ConservationLaw::VariableType::velocity);
    const int imag = law.getIndexFor (ConservationLaw::VariableType::magnetic);
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
            const double vx = 0.05 * std::sin (2 * M_PI * y);
            const double vy = 0.05 * std::cos (2 * M_PI * x);

            for (int k = 0; k < ng; ++k)
            {
                P (i, j, k,           ivel + 0) =  vx;
                P (i, j, k + nk + ng, ivel + 0) = -vx;
                P (i, j, k,           ivel + 1) =  vy;
                P (i, j, k + nk + ng, ivel + 1) = -vy;

                P (i, j, k,           imag + 0) = 0.0;
                P (i, j, k + nk + ng, imag + 0) = 0.0;
                P (i, j, k,           imag + 1) = 0.0;
                P (i, j, k + nk + ng, imag + 1) = 0.0;
            }
        }
    }
}

void DrivenMHDBoundary::applyToCellCenteredB (Cow::Array& B, int numGuard) const
{
    auto periodic = PeriodicBoundaryCondition();
    auto outflow = OutflowBoundaryCondition();

    periodic.applyToAxis (B, numGuard, 0);
    periodic.applyToAxis (B, numGuard, 1);
    outflow.applyToAxis (B, numGuard, 2);
}

void DrivenMHDBoundary::applyToGodunovFluxes (Cow::Array& F, int numGuard, int axis) const
{
    auto periodic = PeriodicBoundaryCondition();
    auto outflow = OutflowBoundaryCondition();

    switch (axis)
    {
        case 0: periodic.applyToAxis (F, numGuard, 0);
        case 1: periodic.applyToAxis (F, numGuard, 1);
        case 2: outflow.applyToAxis (F, numGuard, 2);
    }
}
