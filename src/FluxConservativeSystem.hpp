#ifndef FluxConservativeSystem_hpp
#define FluxConservativeSystem_hpp

#include <functional>
#include "Array.hpp"
#include "Mara.hpp"
#include "ConstrainedTransport.hpp"




class FluxConservativeSystem
{
public:
    /**
    An exception class that holds a list of zones that have failed in a given
    time step.
    */
    class SolverFailure : public std::exception
    {
    public:
        const char* what() const noexcept override;
        std::vector<ConservationLaw::StateFailure> failedStates;
    };

    /**
    Constructor, stores shared pointers to needed algorithms in the given
    setup struct.
    */
    FluxConservativeSystem (SimulationSetup setup);

    /**
    Return a reference to the updateable region in the primitive variable
    array. All fields are returned, but if fieldIndex is non-negative then
    only that component of the primitive data is returned.
    */
    Cow::Array::Reference getPrimitive (int fieldIndex=-1);

    /**
    Return three contiguous components of the primitive data, starting at the
    given index.
    */
    Cow::Array::Reference getPrimitiveVector (int fieldIndex);

    /**
    Return an array indicating the health of corresponding zones. 0 means no
    errors occured there.
    */
    Cow::Array::Reference getZoneHealth();

    /**
    Return the volume integral (not volume average) of named diagnostics
    defined by the ConservationLaw that's in use. This operation uses the
    most recent values of the primitive quantities, as would be retured by
    getPrimitive(). Guard zones are not included in the integral, so the
    result may be added dumbly to that of other grid patches.
    */
    std::vector<double> volumeIntegratedDiagnostics();

    /**
    Set initial data on the primitive variable array.
    */
    void setInitialData (InitialDataFunction F, InitialDataFunction A);

    /**
    Assign a block of primitive data. The size of the array must be the number
    of cells in the mesh, with the size of axis 3 being the number of
    primitive variable fields. This is the same shape as the array returned by
    getPrimitive() when fieldIndex=-1.
    */
    void assignPrimitive (Cow::Array primitiveData);

    /**
    Get the shortest time a wave to cross a cell.
    */
    double getCourantTimestep();

    /**
    Advance the solution in one high level method.
    */
    void advance (double dt);

private:
    /* @internal */
    void computeIntercellFluxes();
    /* @internal */
    void intercellFluxSweep (int axis);
    /* @internal */
    void computeTimeDerivative();
    /* @internal */
    void applyBoundaryCondition();
    /* @internal */
    void cacheConserved();
    /* @internal */
    void updateConserved (double dt);
    /* @internal */
    void averageRungeKutta (double b);
    /* @internal */
    void recoverPrimitive();
    /* @internal */
    void takeRungeKuttaSubstep (double dt, double b);
    /* @internal */
    void uploadFieldsToCT();
    /* @internal */
    UniformCartesianCT* getCT();

    unsigned int numConserved;
    unsigned int stencilSize;
    unsigned int rungeKuttaOrder;
    bool disableCT;

    Cow::Shape domainShape;
    Cow::Region updateableRegion;
    Cow::Region magneticIndices;

    Cow::Array U; // Conserved quantities
    Cow::Array P; // Primitive quantities
    Cow::Array L; // Time derivative of U (to save memory this could be removed)
    Cow::Array U0; // Conserved quantities (cached for RK update)
    Cow::Array F1; // Flux along axis 1
    Cow::Array F2; // Flux along axis 2
    Cow::Array F3; // Flux along axis 3
    Cow::Array zoneHealth; // Nonzero values indicate where failures are occuring

    std::shared_ptr<MeshGeometry> meshGeometry;
    std::shared_ptr<ConservationLaw> conservationLaw;
    std::shared_ptr<IntercellFluxScheme> intercellFluxScheme;
    std::shared_ptr<ConstrainedTransport> constrainedTransport;
    std::shared_ptr<BoundaryCondition> boundaryCondition;
};

#endif
