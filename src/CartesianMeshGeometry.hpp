#include <cassert>
#include "Mara.hpp"




class CartesianMeshGeometry : public MeshGeometry
{
public:
    CartesianMeshGeometry();
    Cow::Shape cellsShape() const override;
    unsigned long totalCellsInMesh() const override;
    Coordinate coordinateAtIndex (double i, double j, double k) const override;
    double faceArea (int i, int j, int k, int axis) const override;
    double cellLength (int i, int j, int k, int axis) const override;
    double cellVolume (int i, int j, int k) const override;
    Cow::Array getPointCoordinates (int axis) const;

    Cow::Shape shape;
    Coordinate lower;
    Coordinate upper;
};
