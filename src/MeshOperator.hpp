#ifndef MeshOperator_hpp
#define MeshOperator_hpp

#include "Array.hpp"
#include "Mara.hpp"




class MeshOperator
{
public:
    using Array = Cow::Array;
    using Shape = Cow::Shape;

    MeshOperator();

    void setMeshGeometry (std::shared_ptr<MeshGeometry>);

    /**
    Return an array which results by applying F to the coordinates (x, y, z)
    of each location in the mesh. The size of the returned array in the first
    three axes is determined by how many of the given type of locations exist
    in the mesh geometry. The size of the returned array on the 4th axis is
    the number of components returned by F, and the 5th axis has length 1.
    */
    Array generate (InitialDataFunction F, MeshLocation location) const;

    Array divergence (const Array& flux, MeshLocation location) const;

    Array curl (const Array& potential, MeshLocation location) const;

private:

    /* @internal */
    Shape makeFacesShape (int numComponents) const;
    /* @internal */
    Shape makeCellsShape (int numComponents) const;

    std::shared_ptr<MeshGeometry> geometry;
    Array faceAreas;
    Array cellVolumes;
};

#endif
