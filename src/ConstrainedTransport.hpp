#ifndef ConstrainedTransport_hpp
#define ConstrainedTransport_hpp


#include "Mara.hpp"




/**
A class that implements CT operations for a special case of uniform, logically
cartesian mesh.
*/
class UniformCartesianCT : public ConstrainedTransport
{
public:
    using Array = Cow::Array;
    using Shape = Cow::Shape;
    using Region = Cow::Region;

    /**
    Data structure returned by methods which request fluxes.
    */
    struct FluxArrays
    {
        Array F1;
        Array F2;
        Array F3;
    };

    /**
    Assign a mesh geometry. This must be done before CT operations are used.
    The mesh geometry is not required at construction, so that algorithm
    choices may be sorted out before the domain size is known.
    */
    void setMeshGeometry (std::shared_ptr<MeshGeometry>) override;

    /**
    Assign a boundary condition to use for Godunov fluxes and magnetic field
    values. This must be done before CT operations are used.
    */
    void setBoundaryCondition (std::shared_ptr<BoundaryCondition>) override;

    /**
    Compute monopole at the given mesh location (vert or cell) from the
    internal array of magnetic field vectors.
    */
    Array computeMonopole (MeshLocation location) const override;

    /**
    Compute the conduction curren (curl B) at the given mesh location from the
    internal array of magnetic field vectors.
    */
    Array computeCurrent (MeshLocation location) const override;

    /**
    Assign Godunov fluxes to the cell faces. Each of the arrays should have 3
    components, although longitudinal magnetic flux, e.g. F1(B1) should always
    be zero.
    */
    void assignGodunovFluxes (Array newF1, Array newF2, Array newF3);

    /**
    Use the given callback function (which returns [Ax, Ay, Az] to assign
    vector potential values at the given mesh location. If location == face,
    then the internal Godunov fluxes are overwritten with the two components
    of A tangent to each face, through the relation Fi(Bj) = -Ek. For this
    operation, vector potential is synonymous with electric field.
    */
    void assignVectorPotential (InitialDataFunction A, MeshLocation location);

    /**
    Assign magnetic field values to cell centers.
    */
    void assignCellCenteredB (Array newB) override;

    /**
    Assign magnetic flux values to faces, formally the dimensions are B * dA,
    but for this class all faces have unit area.
    */
    void assignFaceCenteredH (Array newH);

    /**
    Assign EMF's to cell edges, formally with dimensions E * dx. In general,
    EMF data is likely to be computed internally by CT operations. However,
    assigning to the electric field directly may be useful in generating
    divergenceless magnetic vector fields.
    */
    void assignEdgeCenteredE (Array newE);

    /**
    This operation returns Godunov fluxes based on the Finite Volume Flux-CT
    scheme of Toth (2000) Section 4.4, which may be used to advance the cell-
    centered magnetic field.
    */
    FluxArrays computeGodunovFluxesFieldCT();

    /**
    Retrieve the Godunov fluxes on faces.
    */
    FluxArrays getGodunovFluxes();

private:
    Array computeMonopoleVert() const;
    Array computeMonopoleCell() const;
    void setFaceBC();
    void setCellBC();
    std::shared_ptr<MeshGeometry> meshGeometry;
    std::shared_ptr<BoundaryCondition> boundaryCondition;

    Region updateableRegionF1;
    Region updateableRegionF2;
    Region updateableRegionF3;
    Region updateableRegionB;

    Array F1; // Godunov flux along each axis, stored on faces
    Array F2;
    Array F3;
    Array B; // Magnetic field at cell centers
    Array H; // Magnetic flux through faces
    Array E; // EMF on cell edges (may also be used as vector potential)
};

#endif
