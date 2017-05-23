#include <iostream> // DEBUG
#include <cassert>
#include <cmath>
#include "ConstrainedTransport.hpp"

using namespace Cow;




void UniformCartesianCT::setDomainShape (Cow::Shape shape)
{
    domainShape = shape;
}

Array UniformCartesianCT::computeMonopole (MeshLocation location) const
{
    assert (location == MeshLocation::vert);
    auto shape = domainShape;
    shape[0] += 1;
    shape[1] += 1;
    shape[2] += 1;

    auto divB = Array (shape);

    for (int i = 0; i < B.size(0) - 1; ++i)
    {
        for (int j = 0; j < B.size(1) - 1; ++j)
        {
            double Bi1 = 0.5 * (B(i + 1, j, 0, 0) + B(i + 1, j + 1, 0, 0)); // average Bx_{i+1} in the j-direction
            double Bi0 = 0.5 * (B(i + 0, j, 0, 0) + B(i + 0, j + 1, 0, 0)); // average Bx_{i+0} in the j-direction
            double Bj1 = 0.5 * (B(i, j + 1, 0, 1) + B(i + 1, j + 1, 0, 1)); // average By_{j+1} in the i-direction
            double Bj0 = 0.5 * (B(i, j + 0, 0, 1) + B(i + 1, j + 0, 0, 1)); // average By_{j+0} in the i-direction

            divB (i + 1, j + 1, 0) = (Bi1 - Bi0) + (Bj1 - Bj0);
            divB (i + 1, j + 1, 1) = (Bi1 - Bi0) + (Bj1 - Bj0);
        }
    }

    return divB;
}

void UniformCartesianCT::assignGodunovFluxes (Array newF1, Array newF2, Array newF3)
{
    Shape shapeF1 = domainShape;
    Shape shapeF2 = domainShape;
    Shape shapeF3 = domainShape;

    shapeF1[0] += 1;
    shapeF2[1] += 1;
    shapeF3[2] += 1;
    shapeF1[3] = 3;
    shapeF2[3] = 3;
    shapeF3[3] = 3;

    assert (shapeF1 == newF1.shape());
    assert (shapeF2 == newF2.shape());
    assert (shapeF3 == newF3.shape());

    F1 = std::move (newF1);
    F2 = std::move (newF2);
    F3 = std::move (newF3);
}

void UniformCartesianCT::assignCellCenteredB (Array newB)
{
    assert (newB.size(0) == domainShape[0]);
    assert (newB.size(1) == domainShape[1]);
    assert (newB.size(2) == domainShape[2]);
    assert (newB.size(3) == 3);
    B = std::move (newB);
}

void UniformCartesianCT::assignFaceCenteredH (Array newH)
{
    H = std::move (newH);
}

void UniformCartesianCT::assignEdgeCenteredE (Array newE)
{
    E = std::move (newE);
}

void UniformCartesianCT::computeGodunovFluxesFieldCT (Array& ctF1, Array& ctF2, Array& ctF3) const
{
    ctF1 = F1;
    ctF2 = F2;
    ctF3 = F3;
}
