#include "MeshData.hpp"

using namespace Cow;




MeshData::MeshData (Shape baseShape, Shape boundaryShape, int numComponents) : boundaryShape (boundaryShape)
{
	auto cellsShape = Shape3D (baseShape);

    // Note that interior region is relative, so it can be used on face and
    // cell data even though they have different sizes.

    for (int n = 0; n < 3; ++n)
    {
        cellsShape[n] += 2 * boundaryShape[n];
        interior.lower[n] =  boundaryShape[n];
        interior.upper[n] = -boundaryShape[n];
    }

    P = Array (cellsShape.withRank(1).withComponents (numComponents));
    B = Array (cellsShape.withRank(3).withComponents(1).increased(1));
    Z = Array (cellsShape.withRank(1).withComponents(1));

    velocityIndex = -1;
    magneticIndex = -1;
}

void MeshData::assignPrimitive (Array newP)
{
    P[interior] = newP;
}

void MeshData::assignMagneticField (Array newB, MeshLocation location)
{
    if (magneticIndex == -1)
    {
        throw std::logic_error ("Attempt to assign magnetic field data when magnetic index = -1 (call setMgneticIndex)");
    }
    const int bi = magneticIndex;

    switch (location)
    {
        case MeshLocation::cell: P[interior.withRange (3, bi, bi + 3)] = newB; break;
        case MeshLocation::face: B[interior] = newB; break;
        default: throw std::logic_error ("Bad mesh location for magnetic field data");
    }
}

Array::Reference MeshData::getPrimitive (int i)
{
    return P[i == -1 ? interior : interior.withRange (3, i, i + 1)];
}

Array::Reference MeshData::getPrimitiveVector (int i)
{
    return P[interior.withRange (3, i, i + 3)];
}

Array::Reference MeshData::getMagneticField (MeshLocation location)
{
    if (magneticIndex == -1)
    {
        throw std::logic_error ("Attempt to retrieve non-existent magnetic field data");
    }

    switch (location)
    {
        case MeshLocation::cell: return getPrimitiveVector (magneticIndex);
        case MeshLocation::face: return B[interior];
        default: throw std::logic_error ("Bad mesh location for magnetic field data");
    }
}

Array::Reference MeshData::getZoneHealth()
{
    return Z;
}

Shape3D MeshData::getBoundaryShape() const
{
    return boundaryShape;
}

void MeshData::applyBoundaryCondition (BoundaryCondition& bc)
{
    for (int axis = 0; axis < 3; ++axis)
    {
        if (boundaryShape[axis] > 0)
        {
            bc.apply (P, MeshLocation::cell, MeshBoundary::left , axis, boundaryShape[axis]);
            bc.apply (P, MeshLocation::cell, MeshBoundary::right, axis, boundaryShape[axis]);
        }
    }
}
