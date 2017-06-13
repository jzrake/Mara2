#include "CartesianMeshGeometry.hpp"




// ============================================================================
MeshGeometry::MeshGeometry()
{
    patchIndex[0] = 0;
    patchIndex[1] = 0;
    patchIndex[2] = 0;
    patchIndex[3] = 0;
    patchIndex[4] = 0;
}

std::vector<bool> MeshGeometry::fleshedOutAxes() const
{
    auto shape = cellsShape();
    return { shape[0] > 1, shape[1] > 1, shape[2] > 1};
}

void MeshGeometry::assignPatchIndex (PatchIndex newPatchIndex)
{
    patchIndex = newPatchIndex;
}

MeshGeometry::PatchIndex MeshGeometry::getPatchIndex() const
{
    return patchIndex;
}

MeshGeometry::Coordinate MeshGeometry::coordinateAtIndex (Cow::Index index) const
{
    return coordinateAtIndex (index[0], index[1], index[2]);
}




// ============================================================================
CartesianMeshGeometry::CartesianMeshGeometry()
{
    shape = {{128, 1, 1, 1, 1}};
    lower = {{0.0, 0.0, 0.0}};
    upper = {{1.0, 1.0, 1.0}};
}


void CartesianMeshGeometry::setCellsShape (Cow::Shape S)
{
    shape[0] = S[0];
    shape[1] = S[1];
    shape[2] = S[2];
}

void CartesianMeshGeometry::setLowerUpper (Coordinate L, Coordinate U)
{
    lower = L;
    upper = U;
}


Cow::Shape CartesianMeshGeometry::cellsShape() const
{
    return shape;
}

unsigned long CartesianMeshGeometry::totalCellsInMesh() const
{
    return shape[0] * shape[1] * shape[2];
}

MeshGeometry::Coordinate CartesianMeshGeometry::coordinateAtIndex (double i, double j, double k) const
{
    return Coordinate ({{
        lower[0] + (upper[0] - lower[0]) * (i + 0.5) / shape[0],
        lower[1] + (upper[1] - lower[1]) * (j + 0.5) / shape[1],
        lower[2] + (upper[2] - lower[2]) * (k + 0.5) / shape[2]}});
}

double CartesianMeshGeometry::cellLength (int i, int j, int k, int axis) const
{
    return (upper[axis] - lower[axis]) / shape[axis];
}

double CartesianMeshGeometry::cellVolume (int i, int j, int k) const
{
    const double dx = cellLength (i, j, k, 0);
    const double dy = cellLength (i, j, k, 1);
    const double dz = cellLength (i, j, k, 2);
    return dx * dy * dz;
}

double CartesianMeshGeometry::meshVolume() const
{
    return (upper[0] - lower[0]) * (upper[1] - lower[1]) * (upper[2] - lower[2]);
}

double CartesianMeshGeometry::faceArea (int i, int j, int k, int axis) const
{
    const double dx = cellLength (i, j, k, 0);
    const double dy = cellLength (i, j, k, 1);
    const double dz = cellLength (i, j, k, 2);

    switch (axis)
    {
        case 0: return dy * dz;
        case 1: return dz * dx;
        case 2: return dx * dy;
        default: assert (false);
    }
}

UnitVector CartesianMeshGeometry::faceNormal (int i, int j, int k, int axis) const
{
    switch (axis)
    {
        case 0: return UnitVector::xhat;
        case 1: return UnitVector::yhat;
        case 2: return UnitVector::zhat;
        default: throw std::logic_error ("Axis argument not 0, 1, or 2");
    }
}

double CartesianMeshGeometry::edgeLength (int i, int j, int k, int axis) const
{
    return cellLength (i, j, k, axis);
}

UnitVector CartesianMeshGeometry::edgeVector (int i, int j, int k, int axis) const
{
    switch (axis)
    {
        case 0: return UnitVector::xhat;
        case 1: return UnitVector::yhat;
        case 2: return UnitVector::zhat;
        default: throw std::logic_error ("Axis argument not 0, 1, or 2");
    }
}

Cow::Array CartesianMeshGeometry::getPointCoordinates (int axis) const
{
    assert (0 <= axis && axis < 3);
    auto coords = Cow::Array (shape[axis] + 1);

    for (int n = 0; n < coords.size(); ++n)
    {
        coords[n] = coordinateAtIndex (n - 0.5, n - 0.5, n - 0.5)[axis];
    }
    return coords;
}
