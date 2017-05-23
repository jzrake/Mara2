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
    Constructor takes the number of cells along each axis. Uniform mesh
    spacing is assumed, so no geometrical information is needed.
    */
    UniformCartesianCT (Shape domainShape);

    /**
    Compute monopole level at given mesh location (vert or cell).
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
    void assignCellCenteredB (Array newB);

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
    void computeGodunovFluxesFieldCT (Array& ctF1, Array& ctF2, Array& ctF3) const;

private:
    Shape domainShape;
    Array F1; // Godunov flux along each axis
    Array F2;
    Array F3;
    Array B; // Magnetic field at cell centers
    Array H; // Magnetic flux through faces
    Array E; // EMF on cell edges (may also be used as vector potential)
};

#endif
