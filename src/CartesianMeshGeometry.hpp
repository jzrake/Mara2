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

    Cow::Shape domainShape() const override
    {
        return shape;
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
        return 1.0; // 1D only !!!
    }

    double cellLength (int i, int j, int k, int axis) const override
    {
        return (upper[axis] - lower[axis]) / shape[axis]; // uniform spacing only !!!
    }

    double cellVolume (int i, int j, int k) const override
    {
        return (upper[0] - lower[0]) / shape[0]; // 1D only !!!
    }

    Cow::Shape shape;
    Coordinate lower;
    Coordinate upper;
};
