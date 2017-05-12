#ifndef FluxConservativeSystem_hpp
#define FluxConservativeSystem_hpp

#include <functional>
#include "Array.hpp"
#include "Mara.hpp"




class FluxConservativeSystem
{
public:
    FluxConservativeSystem (SimulationSetup setup);
    Cow::Array::Reference getPrimitive();
    void setInitialData (InitialDataFunction F);
    void computeIntercellFluxes();
    void computeTimeDerivative();
    void applyBoundaryConditions();
    void updateConserved (double dt, double rungeKuttaParameter=0);
    void recoverPrimitive();
    double getCourantTimestep();

private:
    int numDimensions;
    int numConserved;
    int schemeOrder;

    Cow::Shape domainShape;
    Cow::Shape domainShapeAndState;
    Cow::Region updateableRegion;

    Cow::Array U; // Conserved quantities
    Cow::Array P; // Primitive quantities
    Cow::Array L; // Time derivative of U
    Cow::Array F1; // Flux along axis 1
    Cow::Array F2; // Flux along axis 2
    Cow::Array F3; // Flux along axis 3

    std::shared_ptr<MeshGeometry> meshGeometry;
    std::shared_ptr<ConservationLaw> conservationLaw;
    std::shared_ptr<IntercellFluxScheme> intercellFluxScheme;
    std::shared_ptr<BoundaryCondition> boundaryCondition;
};

#endif
