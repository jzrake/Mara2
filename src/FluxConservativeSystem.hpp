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
    Set initial data on the primitive variable array.
    */
    void setInitialData (InitialDataFunction F);

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

    int numConserved;
    int stencilSize;
    int rungeKuttaOrder;

    Cow::Shape domainShape;
    Cow::Region updateableRegion;

    Cow::Array U; // Conserved quantities
    Cow::Array P; // Primitive quantities
    Cow::Array L; // Time derivative of U (to save memory this could be removed)
    Cow::Array U0; // Conserved quantities (cached for RK update)
    Cow::Array F1; // Flux along axis 1
    Cow::Array F2; // Flux along axis 2
    Cow::Array F3; // Flux along axis 3

    std::shared_ptr<MeshGeometry> meshGeometry;
    std::shared_ptr<ConservationLaw> conservationLaw;
    std::shared_ptr<IntercellFluxScheme> intercellFluxScheme;
    std::shared_ptr<ConstrainedTransport> constrainedTransport;
    std::shared_ptr<BoundaryCondition> boundaryCondition;
};

#endif
