#ifndef MeshOperator_hpp
#define MeshOperator_hpp

#include "Array.hpp"
#include "Mara.hpp"




struct GodunovStencil
{
public:
    using Array = Cow::Array;
    using IntercellFluxFunction = std::function<void (GodunovStencil&)>;

    GodunovStencil (int footprint, int numCellQ, int numFaceQ);
    std::vector<Coordinate> cellCoords;
    Array cellData;
    Array faceData;
    Array faceFlux;
    UnitVector faceNormal;
};




class MeshOperator
{
public:
    using Array = Cow::Array;
    using Shape = Cow::Shape;
    using Index = Cow::Index;
    enum class VectorMode { scalars, fluxish, emflike };

    MeshOperator (std::shared_ptr<MeshGeometry> geometry=nullptr);

    void setMeshGeometry (std::shared_ptr<MeshGeometry>);

    /**
    Return an array which results by applying F to the coordinates (x, y, z)
    of each location in the mesh. The size of the returned array in the first
    three axes is determined by how many of the given type of locations exist
    in the mesh geometry. In scalars mode, the size of the returned array on
    the 4th axis is the number of components returned by F. In fluxish mode, F
    must return a 3-component vector, whose projection onto the face normal is
    returned as the sole entry on the 4th axis.
    */
    Array generate (InitialDataFunction F, MeshLocation location,
        VectorMode vectorMode=VectorMode::scalars) const;

    /**
    Compute the divergence of a two-form (magnetic flux densities or flow of
    conserved quantities) on mesh faces. The input array (flux) must have
    size=3 on its 4th axis, and may have any number of components on its 3rd
    axis. One divergence is computed per flux component. The optional start
    index indicates the left-most index of the input data with respect to the
    left of the mesh geometry instance, and must be used if passing in data
    that has guard zones.
    */
    Array divergence (const Array& flux, Index start={}) const;

    /**
    Compute the curl of a one-form on mesh edges. The input array (potential)
    must have size=3 on its 4th axis, and must have one component on its 3rd
    axis. The optional start index indicates the left-most index of the input
    data with respect to the left of the mesh geometry instance, and must be
    used if passing in data that has guard zones.
    */
    Array curl (const Array& potential, Index start={}) const;


    Array godunov (
        GodunovStencil::IntercellFluxFunction Fhat,
        const Array& cellData,
        const Array& faceData,
        Shape footprint,
        Index start={}) const;

private:
    class RichShape;
    std::shared_ptr<MeshGeometry> geometry;
};

#endif
