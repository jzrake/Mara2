#include "CellCenteredFieldCT.hpp"

using namespace Cow;




// ============================================================================
void CellCenteredFieldCT::correctGodunovFluxes (Array& F, int magneticIndex) const
{
    // This function implements the stencil shown in Fig. 3 of Toth (2000), and
    // the formula in Equation (25).

    auto G = F;

    const int B0 = magneticIndex + 0;
    const int B1 = magneticIndex + 1;
    const int B2 = magneticIndex + 2;

    // ------------------------------------------------------------------------
    for (int i = 0; i < F.size(0); ++i)
    for (int j = 0; j < F.size(1); ++j)
    for (int k = 0; k < F.size(2); ++k)
    {
        G (i, j, k, B0, 0) = 0.0;
        G (i, j, k, B1, 1) = 0.0;
        G (i, j, k, B2, 2) = 0.0;
    }

    // ------------------------------------------------------------------------
    if (F.size(0) >= 3)
    {
        for (int i = 1; i < F.size(0) - 0; ++i)
        for (int j = 1; j < F.size(1) - 1; ++j)
        for (int k = 0; k < F.size(2) - 0; ++k)
        {
            G (i, j, k, B1, 0) = 0.125 * (
                + 2 * F(i + 0, j + 0, k + 0, B1, 0) // Fx(By)
                + 1 * F(i + 0, j + 1, k + 0, B1, 0)
                + 1 * F(i + 0, j - 1, k + 0, B1, 0)
                - 1 * F(i + 0, j + 1, k + 0, B0, 1) // Fy(Bx)
                - 1 * F(i - 1, j + 1, k + 0, B0, 1)
                - 1 * F(i + 0, j + 0, k + 0, B0, 1)
                - 1 * F(i - 1, j + 0, k + 0, B0, 1));
        }

        for (int i = 1; i < F.size(0) - 0; ++i)
        for (int j = 0; j < F.size(1) - 0; ++j)
        for (int k = 1; k < F.size(2) - 1; ++k)
        {
            G (i, j, k, B2, 0) = 0.125 * (
                + 2 * F(i + 0, j + 0, k + 0, B2, 0) // Fx(Bz)
                + 1 * F(i + 0, j + 0, k + 1, B2, 0)
                + 1 * F(i + 0, j + 0, k - 1, B2, 0)
                - 1 * F(i + 0, j + 0, k + 1, B0, 2) // Fz(Bx)
                - 1 * F(i - 1, j + 0, k + 1, B0, 2)
                - 1 * F(i + 0, j + 0, k + 0, B0, 2)
                - 1 * F(i - 1, j + 0, k + 0, B0, 2));
        }
    }

    // ------------------------------------------------------------------------
    if (F.size(1) >= 3)
    {
        for (int i = 0; i < F.size(0) - 0; ++i)
        for (int j = 1; j < F.size(1) - 0; ++j)
        for (int k = 1; k < F.size(2) - 1; ++k)
        {
            G (i, j, k, B2, 1) = 0.125 * (
                + 2 * F(i + 0, j + 0, k + 0, B2, 1) // Fy(Bz)
                + 1 * F(i + 0, j + 0, k + 1, B2, 1)
                + 1 * F(i + 0, j + 0, k - 1, B2, 1)
                - 1 * F(i + 0, j + 0, k + 1, B1, 2) // Fz(By)
                - 1 * F(i + 0, j - 1, k + 1, B1, 2)
                - 1 * F(i + 0, j + 0, k + 0, B1, 2)
                - 1 * F(i + 0, j - 1, k + 0, B1, 2));
        }

        for (int i = 1; i < F.size(0) - 1; ++i)
        for (int j = 1; j < F.size(1) - 0; ++j)
        for (int k = 0; k < F.size(2) - 0; ++k)
        {
            G (i, j, k, B0, 1) = 0.125 * (
                + 2 * F(i + 0, j + 0, k + 0, B0, 1) // Fy(Bx)
                + 1 * F(i + 1, j + 0, k + 0, B0, 1)
                + 1 * F(i - 1, j + 0, k + 0, B0, 1)
                - 1 * F(i + 1, j + 0, k + 0, B1, 0) // Fx(By)
                - 1 * F(i + 1, j - 1, k + 0, B1, 0)
                - 1 * F(i + 0, j + 0, k + 0, B1, 0)
                - 1 * F(i + 0, j - 1, k + 0, B1, 0));
        }
    }

    // ------------------------------------------------------------------------
    if (F.size(2) >= 3)
    {
        for (int i = 1; i < F.size(0) - 1; ++i)
        for (int j = 0; j < F.size(1) - 0; ++j)
        for (int k = 1; k < F.size(2) - 0; ++k)
        {
            G (i, j, k, B0, 2) = 0.125 * (
                + 2 * F(i + 0, j + 0, k + 0, B0, 2) // Fz(Bx)
                + 1 * F(i + 1, j + 0, k + 0, B0, 2)
                + 1 * F(i - 1, j + 0, k + 0, B0, 2)
                - 1 * F(i + 1, j + 0, k + 0, B2, 0) // Fx(Bz)
                - 1 * F(i + 1, j + 0, k - 1, B2, 0)
                - 1 * F(i + 0, j + 0, k + 0, B2, 0)
                - 1 * F(i + 0, j + 0, k - 1, B2, 0));
        }

        for (int i = 0; i < F.size(0) - 0; ++i)
        for (int j = 1; j < F.size(1) - 1; ++j)
        for (int k = 1; k < F.size(2) - 0; ++k)
        {
            G (i, j, k, B1, 2) = 0.125 * (
                + 2 * F(i + 0, j + 0, k + 0, B1, 2) // Fz(By)
                + 1 * F(i + 0, j + 1, k + 0, B1, 2)
                + 1 * F(i + 0, j - 1, k + 0, B1, 2)
                - 1 * F(i + 0, j + 1, k + 0, B2, 1) // Fy(Bz)
                - 1 * F(i + 0, j + 1, k - 1, B2, 1)
                - 1 * F(i + 0, j + 0, k + 0, B2, 1)
                - 1 * F(i + 0, j + 0, k - 1, B2, 1));
        }
    }

    // Overwrite old fluxes with the new.
    // ------------------------------------------------------------------------
    F = G;
}

Array CellCenteredFieldCT::vectorPotentialToFluxes (Array A) const
{
    auto F = Array (A.shape3D().withComponents(3).withRank(3));

    F.shape3D().deploy ([&] (int i, int j, int k)
    {
        F (i, j, k, 1, 0) = -A(i, j, k, 2, 0);
        F (i, j, k, 2, 0) = +A(i, j, k, 1, 0);
        F (i, j, k, 2, 1) = -A(i, j, k, 0, 1);
        F (i, j, k, 0, 1) = +A(i, j, k, 2, 1);
        F (i, j, k, 0, 2) = -A(i, j, k, 1, 2);
        F (i, j, k, 1, 2) = +A(i, j, k, 0, 2);
    });
    return F;
}

Array CellCenteredFieldCT::generateGodunovFluxes (const Array& F, int magneticIndex) const
{
    auto G = F;
    correctGodunovFluxes (G, magneticIndex);
    return G;
}

Array CellCenteredFieldCT::monopole (Array B, MeshLocation location) const
{
    auto M = Array (B.shape3D().increased(1).withComponents(1));

    // To deal with flat axes, use variables r, s, t
    // ------------------------------------------------------------------------
    const int r = B.size(0) > 1 ? 1 : 0;
    const int s = B.size(1) > 1 ? 1 : 0;
    const int t = B.size(2) > 1 ? 1 : 0;

    auto S = Shape3D (B.size(0) - r, B.size(1) - s, B.size(2) - t);

    S.deploy ([&] (int i, int j, int k)
    {
        const double B0x00 = B(i + 0, j + 0, k + 0, 0);
        const double B0x01 = B(i + 0, j + 0, k + t, 0);
        const double B0x10 = B(i + 0, j + s, k + 0, 0);
        const double B0x11 = B(i + 0, j + s, k + t, 0);
        const double B1x00 = B(i + r, j + 0, k + 0, 0);
        const double B1x01 = B(i + r, j + 0, k + t, 0);
        const double B1x10 = B(i + r, j + s, k + 0, 0);
        const double B1x11 = B(i + r, j + s, k + t, 0);
        const double B0y00 = B(i + 0, j + 0, k + 0, 1);
        const double B0y01 = B(i + r, j + 0, k + 0, 1);
        const double B0y10 = B(i + 0, j + 0, k + t, 1);
        const double B0y11 = B(i + r, j + 0, k + t, 1);
        const double B1y00 = B(i + 0, j + s, k + 0, 1);
        const double B1y01 = B(i + r, j + s, k + 0, 1);
        const double B1y10 = B(i + 0, j + s, k + t, 1);
        const double B1y11 = B(i + r, j + s, k + t, 1);
        const double B0z00 = B(i + 0, j + 0, k + 0, 2);
        const double B0z01 = B(i + 0, j + s, k + 0, 2);
        const double B0z10 = B(i + r, j + 0, k + 0, 2);
        const double B0z11 = B(i + r, j + s, k + 0, 2);
        const double B1z00 = B(i + 0, j + 0, k + t, 2);
        const double B1z01 = B(i + 0, j + s, k + t, 2);
        const double B1z10 = B(i + r, j + 0, k + t, 2);
        const double B1z11 = B(i + r, j + s, k + t, 2);

        const double Bx0 = B0x00 + B0x01 + B0x10 + B0x11;
        const double Bx1 = B1x00 + B1x01 + B1x10 + B1x11;
        const double By0 = B0y00 + B0y01 + B0y10 + B0y11;
        const double By1 = B1y00 + B1y01 + B1y10 + B1y11;
        const double Bz0 = B0z00 + B0z01 + B0z10 + B0z11;
        const double Bz1 = B1z00 + B1z01 + B1z10 + B1z11;

        M (i + 1, j + 1, k + 1) = 0.25 * ((Bx1 - Bx0) + (By1 - By0) + (Bz1 - Bz0));
    });

    switch (location)
    {
        case MeshLocation::vert: return M;
        case MeshLocation::cell:
        {
            auto Mc = Array (M.shape3D().reduced(1));

            Mc.shape3D().deploy ([&] (int i, int j, int k)
            {
                const double Mi = M (i, j, k) + M (i + 1, j, k);
                const double Mj = M (i, j, k) + M (i, j + 1, k);
                const double Mk = M (i, j, k) + M (i, j, k + 1);
                Mc (i, j, k) = 1. / 6 * (Mi + Mj + Mk);
            });
            return Mc;
        }
        default:
        {
            throw std::logic_error ("Mesh location is not cell or vert");
        }
    }
}

Array CellCenteredFieldCT::current (Array B, MeshLocation location) const
{
    // To deal with flat axes, use variables r, s, t
    // ------------------------------------------------------------------------
    const int r = B.size(0) > 1 ? 1 : 0;
    const int s = B.size(1) > 1 ? 1 : 0;
    const int t = B.size(2) > 1 ? 1 : 0;

    auto S = Shape3D (B.size(0) - 2 * r, B.size(1) - 2 * s, B.size(2) - 2 * t);
    auto J = Array (B.shape());

    S.deploy ([&] (int i, int j, int k)
    {
        const double d0B1 = B (i + r, j, k, 1) - B (i - r, j, k, 1);
        const double d0B2 = B (i + r, j, k, 2) - B (i - r, j, k, 2);
        const double d1B2 = B (i, j + s, k, 2) - B (i, j - s, k, 2);
        const double d1B0 = B (i, j + s, k, 0) - B (i, j - s, k, 0);
        const double d2B0 = B (i, j, k + t, 0) - B (i, j, k - t, 0);
        const double d2B1 = B (i, j, k + t, 1) - B (i, j, k - t, 1);

        J (i + r, j + s, k + t, 0) = d1B2 - d2B1;
        J (i + r, j + s, k + t, 1) = d2B0 - d0B2;
        J (i + r, j + s, k + t, 2) = d0B1 - d1B0;
    });

    return J;
}
