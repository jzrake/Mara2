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

    /**
    This function needs to be called before any of the algorithms are useful.
    The domain shape is not required at construction, so that algorithm
    choices may be sorted out before the domain size is known.
    */
    void setDomainShape (Cow::Shape shape) override;

    /**
    Assign a boundary condition to use for Godunov fluxes and magnetic field
    values.
    */
    void setBoundaryCondition (std::shared_ptr<BoundaryCondition>) override;

    /**
    Compute monopole at the given mesh location (vert or cell) from the given
    array of magnetic field vectors.
    */
    Array computeMonopole (MeshLocation location) const override;

    /**
    Assign Godunov fluxes to the cell faces. Each of the arrays should have 3
    components, although longitudinal magnetic flux, e.g. F1(B1) should always
    be zero.
    */
    void assignGodunovFluxes (Array newF1, Array newF2, Array newF3);

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
    void computeGodunovFluxesFieldCT (Array& ctF1, Array& ctF2, Array& ctF3);

private:
    std::shared_ptr<BoundaryCondition> boundaryCondition;
    Shape domainShape;

    Cow::Region updateableRegionF1;
    Cow::Region updateableRegionF2;
    Cow::Region updateableRegionF3;
    Cow::Region updateableRegionB;

    Array F1; // Godunov flux along each axis
    Array F2;
    Array F3;
    Array B; // Magnetic field at cell centers
    Array H; // Magnetic flux through faces
    Array E; // EMF on cell edges (may also be used as vector potential)
};

#endif
