#include <iostream> // DEBUG
#include <cassert>
#include "ConstrainedTransport.hpp"

using namespace Cow;



UniformCartesianCT::UniformCartesianCT (Shape domainShape) : domainShape (domainShape)
{

}

Array UniformCartesianCT::computeMonopole (MeshLocation location) const
{
    return Array();
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
