#pragma once
#include "Mara.hpp"



class SphericalMeshGeometry : public MeshGeometry
{
public:
    using Shape = Cow::Shape;
    using Index = Cow::Index;
    SphericalMeshGeometry();
    SphericalMeshGeometry (Cow::Shape shape);
    SphericalMeshGeometry (int ni, int nj, int nk);
    void setCellsShape (Cow::Shape S) override;
    void setLowerUpper (Coordinate L, Coordinate U) override;
    Cow::Shape cellsShape() const override;
    Cow::Index indexAtCoordinate (Coordinate x) const override;
    Coordinate coordinateAtIndex (double i, double j, double k) const override;
    unsigned long totalCellsInMesh() const override;
    double cellLength (int i, int j, int k, int axis) const override;
    double cellVolume (int i, int j, int k) const override;
    double meshVolume() const override;
    double faceArea (int i, int j, int k, int axis) const override;
    UnitVector faceNormal (int i, int j, int k, int axis) const override;
    double edgeLength (int i, int j, int k, int axis) const override;
    UnitVector edgeVector (int i, int j, int k, int axis) const override;
    Cow::Array getPointCoordinates (int axis) const;
private:
    void cacheSpacing();
    Coordinate lower;
    Coordinate upper;
    Shape shape;
    std::vector<double> edges[3];
};
