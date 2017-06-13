#include <iostream> // DEBUG
#include <cassert>
#include <cmath>
#include "ConstrainedTransport.hpp"
#include "Stencil.hpp"

using namespace Cow;




void UniformCartesianCT::setMeshGeometry (std::shared_ptr<MeshGeometry> geometry)
{
    meshGeometry = geometry;
    auto domainShape = meshGeometry->cellsShape();

    // The flux arrays we work with will have one more face than cells in the
    // longitudinal direction, and two more faces than cells in the transverse
    // directions.
    auto shapeF1 = domainShape;
    auto shapeF2 = domainShape;
    auto shapeF3 = domainShape;

    shapeF1[0] += 1;
    shapeF1[1] += 2;
    shapeF1[2] += 2;
    shapeF1[3] = 3;
    shapeF2[0] += 2;
    shapeF2[1] += 1;
    shapeF2[2] += 2;
    shapeF2[3] = 3;
    shapeF3[0] += 2;
    shapeF3[1] += 2;
    shapeF3[2] += 1;
    shapeF3[3] = 3;

    F1 = Cow::Array (shapeF1);
    F2 = Cow::Array (shapeF2);
    F3 = Cow::Array (shapeF3);

    updateableRegionF1 = Region();
    updateableRegionF2 = Region();
    updateableRegionF3 = Region();
    updateableRegionF1.lower[1] =  1;
    updateableRegionF1.upper[1] = -1;
    updateableRegionF1.lower[2] =  1;
    updateableRegionF1.upper[2] = -1;
    updateableRegionF2.lower[2] =  1;
    updateableRegionF2.upper[2] = -1;
    updateableRegionF2.lower[0] =  1;
    updateableRegionF2.upper[0] = -1;
    updateableRegionF3.lower[0] =  1;
    updateableRegionF3.upper[0] = -1;
    updateableRegionF3.lower[1] =  1;
    updateableRegionF3.upper[1] = -1;


    // Cell-centered magnetic field has one guard zone on all dimensions
    auto Bshape = domainShape;
    Bshape[0] += 2;
    Bshape[1] += 2;
    Bshape[2] += 2;
    Bshape[3] = 3;

    B = Array (Bshape);

    updateableRegionB = Region();
    updateableRegionB.lower[0] =  1;
    updateableRegionB.upper[0] = -1;
    updateableRegionB.lower[1] =  1;
    updateableRegionB.upper[1] = -1;
    updateableRegionB.lower[2] =  1;
    updateableRegionB.upper[2] = -1;
}

void UniformCartesianCT::setBoundaryCondition (std::shared_ptr<BoundaryCondition> newBC)
{
    boundaryCondition = newBC;
}

Array UniformCartesianCT::computeMonopole (MeshLocation location) const
{
    switch (location)
    {
        case MeshLocation::vert: return computeMonopoleVert();
        case MeshLocation::cell: return computeMonopoleCell();
        default: throw std::logic_error ("Mesh location for monopole calculation not cell or vert");
    }
}

Array UniformCartesianCT::computeCurrent (MeshLocation location) const
{
    if (location != MeshLocation::cell)
    {
        throw std::logic_error ("Mesh location for current calculation not cell");
    }

    auto curl = [] (const Array& B, Array& J)
    {
        const double d1B2 = B (1, 0, 0, 2) - B (-1, 0, 0, 2);
        const double d2B1 = B (0, 1, 0, 1) - B ( 0,-1, 0, 1);
        const double d2B3 = B (0, 1, 0, 3) - B ( 0,-1, 0, 3);
        const double d3B2 = B (0, 0, 1, 2) - B ( 0, 0,-1, 2);
        const double d3B1 = B (0, 0, 1, 1) - B ( 0, 0,-1, 1);
        const double d1B3 = B (1, 0, 0, 3) - B (-1, 0, 0, 3);

        J[0] = d2B3 - d3B2;
        J[1] = d3B1 - d1B3;
        J[2] = d1B2 - d2B1;
    };

    auto stencil = Stencil();
    stencil.setCodomainRank (3, 1);
    stencil.setFootprintLower (-1, -1, -1);
    stencil.setFootprintUpper (+1, +1, +1);

    return stencil.evaluate (curl, B);
}

void UniformCartesianCT::assignGodunovFluxes (Array newF1, Array newF2, Array newF3)
{
    F1[updateableRegionF1] = newF1;
    F2[updateableRegionF2] = newF2;
    F3[updateableRegionF3] = newF3;
    setFaceBC();
}

void UniformCartesianCT::assignVectorPotential (InitialDataFunction A, MeshLocation location)
{
    assert (location == MeshLocation::face);

    if (A (0, 0, 0).size() != 3)
    {
        throw std::runtime_error ("Vector potential function returned vector of length != 3");
    }

    for (int i = 0; i < F1.size(0); ++i)
    {
        for (int j = 0; j < F1.size(1); ++j)
        {
            for (int k = 0; k < F1.size(2); ++k)
            {
                auto X = meshGeometry->coordinateAtIndex (i - 0.5, j - 1, k - 1);
                auto a = A (X[0], X[1], X[2]);

                F1 (i, j, k, 1) = -a[2];
                F1 (i, j, k, 2) = +a[1];
            }
        }
    }

    for (int i = 0; i < F2.size(0); ++i)
    {
        for (int j = 0; j < F2.size(1); ++j)
        {
            for (int k = 0; k < F2.size(2); ++k)
            {
                auto X = meshGeometry->coordinateAtIndex (i - 1, j - 0.5, k - 1);
                auto a = A (X[0], X[1], X[2]);

                F2 (i, j, k, 2) = -a[0];
                F2 (i, j, k, 0) = +a[2];
            }
        }
    }

    for (int i = 0; i < F3.size(0); ++i)
    {
        for (int j = 0; j < F3.size(1); ++j)
        {
            for (int k = 0; k < F3.size(2); ++k)
            {
                auto X = meshGeometry->coordinateAtIndex (i - 1, j - 1, k - 0.5);
                auto a = A (X[0], X[1], X[2]);

                F3 (i, j, k, 0) = -a[1];
                F3 (i, j, k, 1) = +a[0];
            }
        }
    }
    setFaceBC();
}

void UniformCartesianCT::assignCellCenteredB (Array newB)
{
    B[updateableRegionB] = newB;
    setCellBC();
}

void UniformCartesianCT::assignFaceCenteredH (Array newH)
{
    H = std::move (newH);
}

void UniformCartesianCT::assignEdgeCenteredE (Array newE)
{
    E = std::move (newE);
}

UniformCartesianCT::FluxArrays UniformCartesianCT::computeGodunovFluxesFieldCT()
{
    // This function implements the stencil shown in Fig. 3 of Toth (2000), and
    // the formula in Equation (25).

    auto ctFluxes = FluxArrays();
    ctFluxes.F1 = Cow::Array (F1[updateableRegionF1].shape());
    ctFluxes.F2 = Cow::Array (F2[updateableRegionF2].shape());
    ctFluxes.F3 = Cow::Array (F3[updateableRegionF3].shape());

    for (int i = 0; i < ctFluxes.F1.size(0); ++i)
    {
        for (int j = 0; j < ctFluxes.F1.size(1); ++j)
        {
            for (int k = 0; k < ctFluxes.F1.size(2); ++k)
            {
                ctFluxes.F1 (i, j, k, 1) = 0.125 * (
                    + 2 * F1(i + 0, j + 1, k + 1, 1) // Fx(By)
                    + 1 * F1(i + 0, j + 2, k + 1, 1)
                    + 1 * F1(i + 0, j + 0, k + 1, 1)
                    - 1 * F2(i + 0, j + 1, k + 1, 0) // Fy(Bx)
                    - 1 * F2(i + 1, j + 1, k + 1, 0)
                    - 1 * F2(i + 0, j + 0, k + 1, 0)
                    - 1 * F2(i + 1, j + 0, k + 1, 0));

                ctFluxes.F1 (i, j, k, 2) = 0.125 * (
                    + 2 * F1(i + 0, j + 1, k + 1, 2) // Fx(Bz)
                    + 1 * F1(i + 0, j + 1, k + 2, 2)
                    + 1 * F1(i + 0, j + 1, k + 0, 2)
                    - 1 * F3(i + 0, j + 1, k + 1, 0) // Fz(Bx)
                    - 1 * F3(i + 1, j + 1, k + 1, 0)
                    - 1 * F3(i + 0, j + 1, k + 0, 0)
                    - 1 * F3(i + 1, j + 1, k + 0, 0));
            }
        }
    }

    for (int i = 0; i < ctFluxes.F2.size(0); ++i)
    {
        for (int j = 0; j < ctFluxes.F2.size(1); ++j)
        {
            for (int k = 0; k < ctFluxes.F2.size(2); ++k)
            {
                ctFluxes.F2 (i, j, k, 2) = 0.125 * (
                    + 2 * F2(i + 1, j + 0, k + 1, 2) // Fy(Bz)
                    + 1 * F2(i + 1, j + 0, k + 2, 2)
                    + 1 * F2(i + 1, j + 0, k + 0, 2)
                    - 1 * F3(i + 1, j + 0, k + 1, 1) // Fz(By)
                    - 1 * F3(i + 1, j + 1, k + 1, 1)
                    - 1 * F3(i + 1, j + 0, k + 0, 1)
                    - 1 * F3(i + 1, j + 1, k + 0, 1));

                ctFluxes.F2 (i, j, k, 0) = 0.125 * (
                    + 2 * F2(i + 1, j + 0, k + 1, 0) // Fy(Bx)
                    + 1 * F2(i + 2, j + 0, k + 1, 0)
                    + 1 * F2(i + 0, j + 0, k + 1, 0)
                    - 1 * F1(i + 1, j + 0, k + 1, 1) // Fx(By)
                    - 1 * F1(i + 1, j + 1, k + 1, 1)
                    - 1 * F1(i + 0, j + 0, k + 1, 1)
                    - 1 * F1(i + 0, j + 1, k + 1, 1));
            }
        }
    }

    for (int i = 0; i < ctFluxes.F3.size(0); ++i)
    {
        for (int j = 0; j < ctFluxes.F3.size(1); ++j)
        {
            for (int k = 0; k < ctFluxes.F3.size(2); ++k)
            {
                ctFluxes.F3 (i, j, k, 0) = 0.125 * (
                    + 2 * F3(i + 1, j + 1, k + 0, 0) // Fz(Bx)
                    + 1 * F3(i + 2, j + 1, k + 0, 0)
                    + 1 * F3(i + 0, j + 1, k + 0, 0)
                    - 1 * F1(i + 1, j + 1, k + 0, 2) // Fx(Bz)
                    - 1 * F1(i + 1, j + 1, k + 1, 2)
                    - 1 * F1(i + 0, j + 1, k + 0, 2)
                    - 1 * F1(i + 0, j + 1, k + 1, 2));

                ctFluxes.F3 (i, j, k, 1) = 0.125 * (
                    + 2 * F3(i + 1, j + 1, k + 0, 1) // Fz(By)
                    + 1 * F3(i + 1, j + 2, k + 0, 1)
                    + 1 * F3(i + 1, j + 0, k + 0, 1)
                    - 1 * F2(i + 1, j + 1, k + 0, 2) // Fy(Bz)
                    - 1 * F2(i + 1, j + 1, k + 1, 2)
                    - 1 * F2(i + 1, j + 0, k + 0, 2)
                    - 1 * F2(i + 1, j + 0, k + 1, 2));
            }
        }
    }
    return ctFluxes;
}

UniformCartesianCT::FluxArrays UniformCartesianCT::getGodunovFluxes()
{
    auto fluxes = FluxArrays();
    fluxes.F1 = F1[updateableRegionF1];
    fluxes.F2 = F2[updateableRegionF2];
    fluxes.F3 = F3[updateableRegionF3];
    return fluxes;
}

Cow::Array UniformCartesianCT::computeMonopoleVert() const
{
    auto Mshape = meshGeometry->cellsShape();
    Mshape[0] += 1;
    Mshape[1] += 1;
    Mshape[2] += 1;

    auto M = Array (Mshape);

    for (int i = 0; i < M.size(0); ++i)
    {
        for (int j = 0; j < M.size(1); ++j)
        {
            for (int k = 0; k < M.size(2); ++k)
            {
                const double B0x00 = B(i + 0, j + 0, k + 0, 0);
                const double B0x01 = B(i + 0, j + 0, k + 1, 0);
                const double B0x10 = B(i + 0, j + 1, k + 0, 0);
                const double B0x11 = B(i + 0, j + 1, k + 1, 0);
                const double B1x00 = B(i + 1, j + 0, k + 0, 0);
                const double B1x01 = B(i + 1, j + 0, k + 1, 0);
                const double B1x10 = B(i + 1, j + 1, k + 0, 0);
                const double B1x11 = B(i + 1, j + 1, k + 1, 0);

                const double B0y00 = B(i + 0, j + 0, k + 0, 1);
                const double B0y01 = B(i + 1, j + 0, k + 0, 1);
                const double B0y10 = B(i + 0, j + 0, k + 1, 1);
                const double B0y11 = B(i + 1, j + 0, k + 1, 1);
                const double B1y00 = B(i + 0, j + 1, k + 0, 1);
                const double B1y01 = B(i + 1, j + 1, k + 0, 1);
                const double B1y10 = B(i + 0, j + 1, k + 1, 1);
                const double B1y11 = B(i + 1, j + 1, k + 1, 1);

                const double B0z00 = B(i + 0, j + 0, k + 0, 2);
                const double B0z01 = B(i + 0, j + 1, k + 0, 2);
                const double B0z10 = B(i + 1, j + 0, k + 0, 2);
                const double B0z11 = B(i + 1, j + 1, k + 0, 2);
                const double B1z00 = B(i + 0, j + 0, k + 1, 2);
                const double B1z01 = B(i + 0, j + 1, k + 1, 2);
                const double B1z10 = B(i + 1, j + 0, k + 1, 2);
                const double B1z11 = B(i + 1, j + 1, k + 1, 2);

                const double Bx0 = B0x00 + B0x01 + B0x10 + B0x11;
                const double Bx1 = B1x00 + B1x01 + B1x10 + B1x11;
                const double By0 = B0y00 + B0y01 + B0y10 + B0y11;
                const double By1 = B1y00 + B1y01 + B1y10 + B1y11;
                const double Bz0 = B0z00 + B0z01 + B0z10 + B0z11;
                const double Bz1 = B1z00 + B1z01 + B1z10 + B1z11;

                M (i, j, k) = 0.25 * ((Bx1 - Bx0) + (By1 - By0) + (Bz1 - Bz0));
            }
        }
    }
    return M;
}

Cow::Array UniformCartesianCT::computeMonopoleCell() const
{
    auto Mv = computeMonopoleVert();
    auto Mc = Array (Mv.size(0) - 1, Mv.size(1) - 1, Mv.size(2) - 1);

    for (int i = 0; i < Mc.size(0); ++i)
    {
        for (int j = 0; j < Mc.size(1); ++j)
        {
            for (int k = 0; k < Mc.size(2); ++k)
            {
                const double Mi = 0.5 * (Mv (i, j, k) + Mv (i + 1, j, k));
                const double Mj = 0.5 * (Mv (i, j, k) + Mv (i, j + 1, k));
                const double Mk = 0.5 * (Mv (i, j, k) + Mv (i, j, k + 1));
                Mc (i, j, k) = 1. / 3 * (Mi + Mj + Mk);
            }
        }
    }
    return Mc;
}

void UniformCartesianCT::setFaceBC()
{
    auto callBC = [&] (Cow::Array& F, int axis)
    {
        boundaryCondition->apply (F,
            BoundaryCondition::MeshLocation::face,
            BoundaryCondition::MeshBoundary::left,
            axis,
            1);
        boundaryCondition->apply (F,
            BoundaryCondition::MeshLocation::face,
            BoundaryCondition::MeshBoundary::right,
            axis,
            1);
    };
    callBC (F1, 1);
    callBC (F1, 2);
    callBC (F2, 2);
    callBC (F2, 0);
    callBC (F3, 0);
    callBC (F3, 1);
}

void UniformCartesianCT::setCellBC()
{
    auto callBC = [&] (Cow::Array& BB, int axis)
    {
        boundaryCondition->apply (BB,
            BoundaryCondition::MeshLocation::cell,
            BoundaryCondition::MeshBoundary::left,
            axis,
            1);
        boundaryCondition->apply (BB,
            BoundaryCondition::MeshLocation::cell,
            BoundaryCondition::MeshBoundary::right,
            axis,
            1);
    };
    callBC (B, 0);
    callBC (B, 1);
    callBC (B, 2);
}
