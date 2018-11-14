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
    Array godunovFlux;
    UnitVector faceNormal;
};




class MeshOperator
{
public:
    using Array = Cow::Array;
    using Index = Cow::Index;
    using Shape = Cow::Shape;
    using Shape3D = Cow::Shape3D;
    using FluxCorrection = std::function<void (const Array&, Array&)>; // args are (P, Fhat)
    enum class VectorMode { scalars, fluxish, emflike };

    MeshOperator (std::shared_ptr<MeshGeometry> geometry=nullptr);

    void setMeshGeometry (std::shared_ptr<MeshGeometry>);

    const MeshGeometry& getMeshGeometry() const;

    /**
    Return an array which results by applying F to the coordinates (x, y, z)
    of each location in the mesh. The size of the returned array in the first
    three axes is determined by how many of the given type of locations exist
    in the mesh geometry. In scalars mode, the size of the returned array on
    the 4th axis is the number of components returned by F. In fluxish mode, F
    must return a 3-component vector, whose projection onto the face normal is
    returned as the sole entry on the 4th axis. An optional boundaryShape
    argument specifies a number of cells to add above the mesh geometry shape
    on each axis.
    */
    Array generate (InitialDataFunction F, MeshLocation location,
        VectorMode vectorMode=VectorMode::scalars, Shape3D boundaryShape=Shape3D()) const;

    /**
    Generate an array of source terms from the given callback function and an
    array of primitive data. The optional start index indicates the left-most
    index of the input data with respect to the left of the mesh geometry
    instance, and must be used if passing in data that has guard zones.
    */
    Array generateSourceTerms (SourceTermsFunction S, const Array& P, Index start={}) const;

    /**
    Return the 1, 2, or 3 dimensional measure of the given mesh locations:
    length of edges, area of faces, or volume of cells.
    */
    Array measure (MeshLocation location) const;

    /**
    Return some measure of the smallest linear dimension in each cell.
    */
    Array linearCellDimension() const;

    /**
    Return an array of the cell centroid coordinates.
    */
    Array cellCentroidCoordinates() const;

    /**
    Compute the divergence of a two-form (magnetic flux densities or flow of
    conserved quantities) on mesh faces. The input array (flux) must have
    size=3 on its 4th axis, and may have any number of components on its 3rd
    axis. One divergence is computed per flux component. The resulting
    divergence is multiplied by the optional factor parameter before being
    returned. The optional start index indicates the left-most index of the
    input data with respect to the left of the mesh geometry instance, and
    must be used if passing in data that has guard zones.
    */
    Array divergence (const Array& flux, double factor=1.0, Index start={}) const;

    /**
    Compute the curl of a one-form on mesh edges. The input array (potential)
    must have size=3 on its 4th axis, and must have one component on its 3rd
    axis. The optional start index indicates the left-most index of the input
    data with respect to the left of the mesh geometry instance, and must be
    used if passing in data that has guard zones.
    */
    Array curl (const Array& potential, Index start={}) const;

    /**
    Return an array of fluxes on cell faces from a dimensionally split Godunov
    operator, Fhat. The footprint argument determines the stencil size, and
    the stencil shape is like a plus sign. Footprint must be an even number,
    as it corresponds to the total number of cells surrounding the face.
    Different axes may have different footprint sizes, but Fhat is used to
    compute the Godunov flux on each axis. The start parameter indicates the
    lower index of the input data, relative to the index origin of the
    geometry instance. If the callback function fluxFunction is provided, that
    function will be called on the flux array before returning it. This can be
    useful, for example in field CT methods which average fluxes of magnetic
    field to keep cell-centered div.B = 0.
    */
    Array godunov (
        GodunovStencil::IntercellFluxFunction Fhat,
        const Array& cellData,
        const Array& faceData,
        Shape footprint,
        Index start={},
        FluxCorrection fluxCorrection=nullptr) const;

    /**
    Assign weight densities to the grid cells based on a coordinate. The
    vectors indexes will contain the mesh indexes over which the input is
    distributed. The weights vector contains the mass in each cell divided by
    the cell volume, such that the sum of volume times assigned density is
    equal to 1. If those vectors already have same size as the number of cells
    receiving weight, then for performance reasons they will not be re-
    allocated.
    */
    void weight (Coordinate point, std::vector<Index>& indexes, std::vector<double>& weights) const;

private:
    std::shared_ptr<MeshGeometry> geometry;
};

#endif
