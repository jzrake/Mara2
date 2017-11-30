#include <cassert>
#include <cmath>
#include "SphericalMeshGeometry.hpp"




// ============================================================================
SphericalMeshGeometry::SphericalMeshGeometry()
{
    shape = {{ 128, 1, 1, 1, 1 }};
    lower = {{ 1.0, 0.0, 0.0 }};
    upper = {{ 10.0, M_PI, 2 * M_PI }};
    cacheSpacing();
}

SphericalMeshGeometry::SphericalMeshGeometry(Cow::Shape S)
{
    shape = {{ S[0], S[1], S[2], 1, 1 }};
    lower = {{ 1.0, 0.0, 0.0 }};
    upper = {{ 10.0, M_PI, 2 * M_PI }};
    cacheSpacing();
}

SphericalMeshGeometry::SphericalMeshGeometry (int ni, int nj, int nk)
{
    shape = {{ ni, nj, nk, 1, 1 }};
    lower = {{ 1.0, 0.0, 0.0 }};
    upper = {{ 10.0, M_PI, 2 * M_PI }};
    cacheSpacing();
}

void SphericalMeshGeometry::setCellsShape (Cow::Shape S)
{
    shape[0] = S[0];
    shape[1] = S[1];
    shape[2] = S[2];
    cacheSpacing();
}

void SphericalMeshGeometry::setLowerUpper (Coordinate L, Coordinate U)
{
    lower = L;
    upper = U;
    cacheSpacing();
}

Cow::Shape SphericalMeshGeometry::cellsShape() const
{
    return shape;
}

Cow::Index SphericalMeshGeometry::indexAtCoordinate (Coordinate x) const
{
    throw std::logic_error ("SphericalMeshGeometry::indexAtCoordinate not implemented");
    return Cow::Index();
}

Coordinate SphericalMeshGeometry::coordinateAtIndex (double i, double j, double k) const
{
    const double r0 = edges[0].at (int(i) + 0);
    const double r1 = edges[0].at (int(i) + 1);
    const double q0 = edges[1].at (int(j) + 0);
    const double q1 = edges[1].at (int(j) + 1);
    const double p0 = edges[2].at (int(k) + 0);
    const double p1 = edges[2].at (int(k) + 1);
    const double r = std::sqrt (r0 * r1) * std::pow (r1 / r0, i);
    const double q = (0.5 - j) * q0 + (0.5 + j) * q1;
    const double p = (0.5 - k) * p0 + (0.5 + k) * p1;
    return Coordinate ({{ r, q, p }});
}

unsigned long SphericalMeshGeometry::totalCellsInMesh() const
{
    return shape[0] * shape[1] * shape[2];
}

double SphericalMeshGeometry::cellLength (int i, int j, int k, int axis) const
{
    const double ri = std::sqrt (edges[0][i] * edges[0][i]);
    const double qi = 0.5 *     (edges[1][j] + edges[1][j]);

    switch (axis)
    {
        case 0: return 1. * (edges[0][i + 1] - edges[0][i]);
        case 1: return ri * (edges[1][j + 1] - edges[1][j]);
        case 2: return ri * (edges[2][k + 1] - edges[2][k]) * std::sin (qi);
        default: throw std::logic_error ("Axis argument not 0, 1, or 2");
    }
}

double SphericalMeshGeometry::cellVolume (int i, int j, int k) const
{
    const double r0 = edges[0].at (i + 0);
    const double r1 = edges[0].at (i + 1);
    const double q0 = edges[1].at (j + 0);
    const double q1 = edges[1].at (j + 1);
    const double p0 = edges[2].at (k + 0);
    const double p1 = edges[2].at (k + 1);
    return -1. / 3 * (r1 * r1 * r1 - r0 * r0 * r0) * (std::cos (q1) - std::cos (q0)) * (p1 - p0);
}

double SphericalMeshGeometry::meshVolume() const
{
    const double r0 = edges[0].front();
    const double r1 = edges[0].back();
    const double q0 = edges[1].front();
    const double q1 = edges[1].back();
    const double p0 = edges[2].front();
    const double p1 = edges[2].back();
    return -1. / 3 * (r1 * r1 * r1 - r0 * r0 * r0) * (std::cos (q1) - std::cos (q0)) * (p1 - p0);
}

double SphericalMeshGeometry::faceArea (int i, int j, int k, int axis) const
{
    const double r0 = edges[0].at (i + 0);
    const double r1 = edges[0].at (i + 1);
    const double q0 = edges[1].at (j + 0);
    const double q1 = edges[1].at (j + 1);
    const double p0 = edges[2].at (k + 0);
    const double p1 = edges[2].at (k + 1);

    switch (axis)
    {
        case 0: return -r0 * r0 * (p1 - p0) * (std::cos (q1) - std::cos (q0));
        case 1: return +1. / 3 * (r1 * r1 * r1 - r0 * r0 * r0) * std::sin (q0) * (p1 - p0);
        case 2: return -1. / 3 * (r1 * r1 * r1 - r0 * r0 * r0) * (std::cos (q1) - std::cos (q0));
        default: throw std::logic_error ("Invalid axis");
    }
}

UnitVector SphericalMeshGeometry::faceNormal (int i, int j, int k, int axis) const
{
    switch (axis)
    {
        case 0: return UnitVector::xhat;
        case 1: return UnitVector::yhat;
        case 2: return UnitVector::zhat;
        default: throw std::logic_error ("Axis argument not 0, 1, or 2");
    }
}

double SphericalMeshGeometry::edgeLength (int i, int j, int k, int axis) const
{
    const double r0 = edges[0][i];
    const double q0 = edges[1][j];

    switch (axis)
    {
        case 0: return 1. * (edges[0][i + 1] - edges[0][i]);
        case 1: return r0 * (edges[1][j + 1] - edges[1][j]);
        case 2: return r0 * (edges[2][k + 1] - edges[2][k]) * std::sin (q0);
        default: throw std::logic_error ("Axis argument not 0, 1, or 2");
    }
}

UnitVector SphericalMeshGeometry::edgeVector (int i, int j, int k, int axis) const
{
    switch (axis)
    {
        case 0: return UnitVector::xhat;
        case 1: return UnitVector::yhat;
        case 2: return UnitVector::zhat;
        default: throw std::logic_error ("Axis argument not 0, 1, or 2");
    }
}

Cow::Array SphericalMeshGeometry::getPointCoordinates (int axis) const
{
    if (axis < 0 || 2 < axis) throw std::logic_error ("Invalid axis");
    auto coords = Cow::Array (shape[axis] + 1);

    for (int n = 0; n < coords.size(); ++n)
    {
        coords[n] = coordinateAtIndex (n - 0.5, n - 0.5, n - 0.5)[axis];
    }
    return coords;
}

void SphericalMeshGeometry::cacheSpacing()
{
    edges[0] = std::vector<double> (shape[0] + 1);
    edges[1] = std::vector<double> (shape[1] + 1);
    edges[2] = std::vector<double> (shape[2] + 1);

    for (int i = 0; i < edges[0].size(); ++i)
    {
        edges[0][i] = lower[0] * std::pow (upper[0] / lower[0], double(i) / edges[0].size());
    }
    for (int j = 0; j < edges[0].size(); ++j)
    {
        edges[1][j] = lower[1] + (upper[1] - lower[1]) * double(j) / edges[1].size();
    }
    for (int k = 0; k < edges[0].size(); ++k)
    {
        edges[2][k] = lower[2] + (upper[2] - lower[2]) * double(k) / edges[2].size();
    }
}
