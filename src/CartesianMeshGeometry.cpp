#include <cassert>
#include "CartesianMeshGeometry.hpp"




// ============================================================================
CartesianMeshGeometry::CartesianMeshGeometry()
{
    shape = {{128, 1, 1, 1, 1}};
    lower = {{0.0, 0.0, 0.0}};
    upper = {{1.0, 1.0, 1.0}};
    cacheSpacing();
}

CartesianMeshGeometry::CartesianMeshGeometry(Cow::Shape S)
{
    shape = {{S[0], S[1], S[2], 1, 1}};
    lower = {{0.0, 0.0, 0.0}};
    upper = {{1.0, 1.0, 1.0}};
    cacheSpacing();
}

CartesianMeshGeometry::CartesianMeshGeometry (int ni, int nj, int nk)
{
    shape = {{ni, nj, nk, 1, 1}};
    lower = {{0.0, 0.0, 0.0}};
    upper = {{1.0, 1.0, 1.0}};
    cacheSpacing();
}

void CartesianMeshGeometry::setCellsShape (Cow::Shape S)
{
    shape[0] = S[0];
    shape[1] = S[1];
    shape[2] = S[2];
    cacheSpacing();
}

void CartesianMeshGeometry::setLowerUpper (Coordinate L, Coordinate U)
{
    lower = L;
    upper = U;
    cacheSpacing();
}

Cow::Shape CartesianMeshGeometry::cellsShape() const
{
    return shape;
}

Cow::Index CartesianMeshGeometry::indexAtCoordinate (Coordinate x) const
{
    return Cow::Index ({{
        int((x[0] - lower[0]) / dx[0]),
        int((x[1] - lower[1]) / dx[1]),
        int((x[2] - lower[2]) / dx[2]),
        0, 0 }});
}

Coordinate CartesianMeshGeometry::coordinateAtIndex (double i, double j, double k) const
{
    return Coordinate ({{
        lower[0] + dx[0] * (i + 0.5),
        lower[1] + dx[1] * (j + 0.5),
        lower[2] + dx[2] * (k + 0.5)}});
}

unsigned long CartesianMeshGeometry::totalCellsInMesh() const
{
    return shape[0] * shape[1] * shape[2];
}

double CartesianMeshGeometry::cellLength (int i, int j, int k, int axis) const
{
    return dx[axis];
}

double CartesianMeshGeometry::cellVolume (int i, int j, int k) const
{
    return dV;
}

double CartesianMeshGeometry::meshVolume() const
{
    return (upper[0] - lower[0]) * (upper[1] - lower[1]) * (upper[2] - lower[2]);
}

double CartesianMeshGeometry::faceArea (int i, int j, int k, int axis) const
{
    return dA[axis];
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

std::shared_ptr<MeshGeometry> CartesianMeshGeometry::duplicate() const
{
    auto mg = new CartesianMeshGeometry;
    *mg = *this;
    return std::shared_ptr<MeshGeometry> (mg);
}

std::string CartesianMeshGeometry::getType() const
{
    return "cartesian";
}

void CartesianMeshGeometry::cacheSpacing()
{
    for (int n = 0; n < 3; ++n)
    {
        dx[n] = (upper[n] - lower[n]) / shape[n];
    }

    dA[0] = dx[1] * dx[2];
    dA[1] = dx[2] * dx[0];
    dA[2] = dx[0] * dx[1];

    dV = dx[0] * dx[1] * dx[2];
}
