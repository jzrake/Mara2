#include <cassert>
#include "Mara.hpp"




// ============================================================================
class CartesianMeshGeometry : public MeshGeometry
{
public:
    CartesianMeshGeometry()
    {
        shape = {{128, 1, 1, 1, 1}};
        lower = {{0.0, 0.0, 0.0}};
        upper = {{1.0, 1.0, 1.0}};
    }

    Cow::Shape cellsShape() const override
    {
        return shape;
    }

    unsigned long totalCellsInMesh() const override
    {
        return shape[0] * shape[1] * shape[2];
    }

    Coordinate coordinateAtIndex (double i, double j, double k) const override
    {
        return Coordinate ({{
            lower[0] + (upper[0] - lower[0]) * (i + 0.5) / shape[0],
            lower[1] + (upper[1] - lower[1]) * (j + 0.5) / shape[1],
            lower[2] + (upper[2] - lower[2]) * (k + 0.5) / shape[2]}});
    }

    double faceArea (int i, int j, int k, int axis) const override
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

    double cellLength (int i, int j, int k, int axis) const override
    {
        return (upper[axis] - lower[axis]) / shape[axis]; // uniform spacing only !!!
    }

    double cellVolume (int i, int j, int k) const override
    {
        const double dx = cellLength (i, j, k, 0);
        const double dy = cellLength (i, j, k, 1);
        const double dz = cellLength (i, j, k, 2);
        return dx * dy * dz;
    }

    Cow::Shape shape;
    Coordinate lower;
    Coordinate upper;
};
