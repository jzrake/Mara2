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

void MeshData::assignMagneticField (Array newB, MeshLocation location, int flags)
{
    if (magneticIndex == -1)
    {
        throw std::logic_error ("Attempt to assign magnetic field data when magnetic index = -1 (call setMgneticIndex)");
    }
    const int bi = magneticIndex;

    switch (location)
    {
        case MeshLocation::cell: P[getRegion (flags).withRange (3, bi, bi + 3)] = newB; break;
        case MeshLocation::face: B[getRegion (flags)] = newB; break;
        default: throw std::logic_error ("Bad mesh location for magnetic field data");
    }
}

void MeshData::allocateDiagnostics (std::vector<std::string> diagnosticFieldNamesToUse)
{
    diagnosticFieldNames = diagnosticFieldNamesToUse;
    D = Array (P.shape3D().withComponents (diagnosticFieldNames.size()));
}

void MeshData::assignDiagnostic (Array newD, int index, int flags)
{
    D[getRegion (flags).withRange (3, index, index + 1)] = newD;
}

Array::Reference MeshData::getPrimitive (int i, int flags)
{
    auto R = getRegion (flags);
    return P[i == -1 ? R : R.withRange (3, i, i + 1)];
}

Array::Reference MeshData::getPrimitiveVector (int i, int flags)
{
    return P[getRegion (flags).withRange (3, i, i + 3)];
}

Array::Reference MeshData::getMagneticField (MeshLocation location, int flags)
{
    if (magneticIndex == -1)
    {
        throw std::logic_error ("Attempt to retrieve non-existent magnetic field data");
    }

    switch (location)
    {
        case MeshLocation::cell: return getPrimitiveVector (magneticIndex, flags);
        case MeshLocation::face: return B[getRegion (flags)];
        default: throw std::logic_error ("Bad mesh location for magnetic field data");
    }
}

Array::Reference MeshData::getZoneHealth (int flags)
{
    return Z;
}

Array::Reference MeshData::getDiagnostic (int index, int flags)
{
    return D[getRegion (flags).withRange (3, index, index + 1)];
}

int MeshData::getNumDiagnostics() const
{
    return diagnosticFieldNames.size();
}

std::string MeshData::getDiagnosticName (int index) const
{
    return diagnosticFieldNames[index];
}

Shape3D MeshData::getBoundaryShape() const
{
    return boundaryShape;
}

void MeshData::applyBoundaryCondition (BoundaryCondition& bc)
{
    bc.applySimple (P, boundaryShape);
}

Region MeshData::getRegion (int flags) const
{
    return flags & includeGuard ? Region() : interior;
}
