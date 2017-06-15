#include <cassert>
#include "Mara.hpp"




class CartesianMeshGeometry : public MeshGeometry
{
public:
    CartesianMeshGeometry();

    CartesianMeshGeometry(int ni, int nj, int nk);
    
    void setCellsShape (Cow::Shape S) override;
    
    void setLowerUpper (Coordinate L, Coordinate U) override;
    
    Cow::Shape cellsShape() const override;
    
    unsigned long totalCellsInMesh() const override;
    
    Coordinate coordinateAtIndex (double i, double j, double k) const override;
    
    double cellLength (int i, int j, int k, int axis) const override;
    
    double cellVolume (int i, int j, int k) const override;
    
    double meshVolume() const override;
    
    double faceArea (int i, int j, int k, int axis) const override;
    
    UnitVector faceNormal (int i, int j, int k, int axis) const override;

    double edgeLength (int i, int j, int k, int axis) const override;

    UnitVector edgeVector (int i, int j, int k, int axis) const override;

    Cow::Array getPointCoordinates (int axis) const;

// private:
    Cow::Shape shape;
    Coordinate lower;
    Coordinate upper;
};
