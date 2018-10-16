#include <cassert>
#include <cmath>
#include "SphericalMeshGeometry.hpp"




// ============================================================================
SphericalMeshGeometry::SphericalMeshGeometry()
{
    shape = {{ 128, 1, 1, 1, 1 }};
    lower = {{ 1.0, 0.0, 0.0 }};
    upper = {{ 10.0, M_PI, 2 * M_PI }};
}

SphericalMeshGeometry::SphericalMeshGeometry(Cow::Shape S)
{
    shape = {{ S[0], S[1], S[2], 1, 1 }};
    lower = {{ 1.0, 0.0, 0.0 }};
    upper = {{ 10.0, M_PI, 2 * M_PI }};
}

SphericalMeshGeometry::SphericalMeshGeometry (int ni, int nj, int nk)
{
    shape = {{ ni, nj, nk, 1, 1 }};
    lower = {{ 1.0, 0.0, 0.0 }};
    upper = {{ 10.0, M_PI, 2 * M_PI }};
}

void SphericalMeshGeometry::setUseLogarithmicRadialBinning (bool shouldUseLogarithmicRadialBinning)
{
    logarithmic = shouldUseLogarithmicRadialBinning;
}

void SphericalMeshGeometry::setCellsShape (Cow::Shape S)
{
    shape[0] = S[0];
    shape[1] = S[1];
    shape[2] = S[2];
}

void SphericalMeshGeometry::setLowerUpper (Coordinate L, Coordinate U)
{
    lower = L;
    upper = U;
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

    const double r = getEdge (i + 0.5, 0);
    const double q = getEdge (j + 0.5, 1);
    const double p = getEdge (k + 0.5, 2);
    return Coordinate ({{ r, q, p }});
}

unsigned long SphericalMeshGeometry::totalCellsInMesh() const
{
    return shape[0] * shape[1] * shape[2];
}

double SphericalMeshGeometry::cellLength (int i, int j, int k, int axis) const
{
    const double ri = logarithmic
    ? std::sqrt (getEdge (i, 0) * getEdge (i + 1, 0))
    : 0.5     * (getEdge (i, 0) + getEdge (i + 1, 1));

    const double qi = 0.5 * (getEdge (j, 1) + getEdge (j + 1, 1));

    switch (axis)
    {
        case 0: return 1. * (getEdge (i + 1, 0) - getEdge (i, 0));
        case 1: return ri * (getEdge (j + 1, 1) - getEdge (j, 1));
        case 2: return ri * (getEdge (k + 1, 2) - getEdge (k, 2)) * std::sin (qi);
        default: throw std::logic_error ("Axis argument not 0, 1, or 2");
    }
}

double SphericalMeshGeometry::cellVolume (int i, int j, int k) const
{
    const double r0 = getEdge (i + 0, 0);
    const double r1 = getEdge (i + 1, 0);
    const double q0 = getEdge (j + 0, 1);
    const double q1 = getEdge (j + 1, 1);
    const double p0 = getEdge (k + 0, 2);
    const double p1 = getEdge (k + 1, 2);
    return -1. / 3 * (r1 * r1 * r1 - r0 * r0 * r0) * (std::cos (q1) - std::cos (q0)) * (p1 - p0);
}

double SphericalMeshGeometry::meshVolume() const
{
    const double r0 = lower[0];
    const double r1 = upper[0];
    const double q0 = lower[1];
    const double q1 = upper[1];
    const double p0 = lower[2];
    const double p1 = upper[2];
    return -1. / 3 * (r1 * r1 * r1 - r0 * r0 * r0) * (std::cos (q1) - std::cos (q0)) * (p1 - p0);
}

double SphericalMeshGeometry::faceArea (int i, int j, int k, int axis) const
{
    const double r0 = getEdge (i + 0, 0);
    const double r1 = getEdge (i + 1, 0);
    const double q0 = getEdge (j + 0, 1);
    const double q1 = getEdge (j + 1, 1);
    const double p0 = getEdge (k + 0, 2);
    const double p1 = getEdge (k + 1, 2);

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
    const double r0 = getEdge (i, 0);
    const double q0 = getEdge (j, 1);

    switch (axis)
    {
        case 0: return 1. * (getEdge (i + 1, 0) - getEdge (i, 0));
        case 1: return r0 * (getEdge (j + 1, 1) - getEdge (j, 1));
        case 2: return r0 * (getEdge (k + 1, 2) - getEdge (k, 2)) * std::sin (q0);
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

std::shared_ptr<MeshGeometry> SphericalMeshGeometry::duplicate() const
{
    auto mg = new SphericalMeshGeometry;
    *mg = *this;
    return std::shared_ptr<MeshGeometry> (mg);
}

std::string SphericalMeshGeometry::getType() const
{
    return "spherical";
}

double SphericalMeshGeometry::getEdge (double n, int axis) const
{
    switch (axis)
    {
        case 0:
        {
            if (logarithmic)
                return lower[0] * std::pow (upper[0] / lower[0], n / shape[0]);
            else
                return lower[0] + (upper[0] - lower[0]) * n / shape[0];
        }
        case 1: return lower[1] + (upper[1] - lower[1]) * n / shape[1];
        case 2: return lower[2] + (upper[2] - lower[2]) * n / shape[2];
        default: assert (false);
    }
}
