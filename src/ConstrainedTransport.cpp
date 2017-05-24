#include <iostream> // DEBUG
#include <cassert>
#include <cmath>
#include "ConstrainedTransport.hpp"

using namespace Cow;




void UniformCartesianCT::setDomainShape (Cow::Shape shape)
{
    domainShape = shape;


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
    assert (location == MeshLocation::vert);

    auto Mshape = domainShape;
    Mshape[0] += 1;
    Mshape[1] += 1;
    Mshape[2] += 1;

    auto M = Array (Mshape);

    for (int i = 0; i < B.size(0) - 1; ++i)
    {
        for (int j = 0; j < B.size(1) - 1; ++j)
        {
            double Bi1 = 0.5 * (B(i + 1, j, 0, 0) + B(i + 1, j + 1, 0, 0)); // average Bx_{i+1} in the j-direction
            double Bi0 = 0.5 * (B(i + 0, j, 0, 0) + B(i + 0, j + 1, 0, 0)); // average Bx_{i+0} in the j-direction
            double Bj1 = 0.5 * (B(i, j + 1, 0, 1) + B(i + 1, j + 1, 0, 1)); // average By_{j+1} in the i-direction
            double Bj0 = 0.5 * (B(i, j + 0, 0, 1) + B(i + 1, j + 0, 0, 1)); // average By_{j+0} in the i-direction

            M (i, j, 0) = (Bi1 - Bi0) + (Bj1 - Bj0);
            M (i, j, 1) = (Bi1 - Bi0) + (Bj1 - Bj0);
        }
    }
    return M;
}

void UniformCartesianCT::assignGodunovFluxes (Array newF1, Array newF2, Array newF3)
{
    F1[updateableRegionF1] = newF1;
    F2[updateableRegionF2] = newF2;
    F3[updateableRegionF3] = newF3;

    boundaryCondition->applyToGodunovFluxes (F1, 1, 1);
    boundaryCondition->applyToGodunovFluxes (F1, 1, 2);
    boundaryCondition->applyToGodunovFluxes (F2, 1, 2);
    boundaryCondition->applyToGodunovFluxes (F2, 1, 0);
    boundaryCondition->applyToGodunovFluxes (F3, 1, 0);
    boundaryCondition->applyToGodunovFluxes (F3, 1, 1);
}

void UniformCartesianCT::assignCellCenteredB (Array newB)
{
    B[updateableRegionB] = newB;
    boundaryCondition->applyToCellCenteredB (B, 1);
}

void UniformCartesianCT::assignFaceCenteredH (Array newH)
{
    H = std::move (newH);
}

void UniformCartesianCT::assignEdgeCenteredE (Array newE)
{
    E = std::move (newE);
}

void UniformCartesianCT::computeGodunovFluxesFieldCT (Array& ctF1, Array& ctF2, Array& ctF3)
{
    // This function implements the stencil shown in Fig. 3 of Toth (2000), and
    // the formula in Equation (25).

    ctF1 = F1[updateableRegionF1];
    ctF2 = F2[updateableRegionF2];
    ctF3 = F3[updateableRegionF3];

    // The work and output flux arrays for Fx have the same size in the x-axis
    assert (ctF1.size(0) == F1.size(0));

    for (int i = 0; i < ctF1.size(0); ++i)
    {
        for (int j = 0; j < ctF1.size(1); ++j)
        {
            ctF1 (i, j, 0, 1) = 0.125 * (0.0
                + 2 * F1(i + 0, j + 1, 0, 1)
                + 1 * F1(i + 0, j + 2, 0, 1)
                + 1 * F1(i + 0, j + 0, 0, 1)
                - 1 * F2(i + 0, j + 1, 0, 0)
                - 1 * F2(i + 1, j + 1, 0, 0)
                - 1 * F2(i + 0, j + 0, 0, 0)
                - 1 * F2(i + 1, j + 0, 0, 0));
        }
    }

    // The work and output flux arrays for Fy have the same size in the y-axis
    assert (ctF2.size(1) == F2.size(1));

    for (int j = 0; j < ctF2.size(1); ++j)
    {
        for (int i = 0; i < ctF2.size(0); ++i)
        {
            ctF2 (i, j, 0, 0) = 0.125 * (0.0   // Fy_{i, j+1/2}
                + 2 * F2(i + 1, j + 0, 0, 0)   // Fy_{i, j+1/2}
                + 1 * F2(i + 2, j + 0, 0, 0)   // Fy_{i+1, j+1/2}
                + 1 * F2(i + 0, j + 0, 0, 0)   // Fy_{i-1, j+1/2}
                - 1 * F1(i + 1, j + 0, 0, 1)   // Fx_{i+1/2, j}
                - 1 * F1(i + 1, j + 1, 0, 1)   // Fx_{i+1/2, j+1}
                - 1 * F1(i + 0, j + 0, 0, 1)   // Fx_{i-1/2, j}
                - 1 * F1(i + 0, j + 1, 0, 1)); // Fx_{i-1/2, j+1}
        }
    }
}
