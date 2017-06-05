#ifndef SimulationSetup_hpp
#define SimulationSetup_hpp

#include "Mara.hpp"
#include "Variant.hpp"




class SimulationSetup2
{
public:
    SimulationSetup2 (Variant::NamedValues config);

    std::shared_ptr<MaraDriver> maraDriver;
    std::shared_ptr<MeshDecomposition> meshDecomposition;
    std::shared_ptr<MeshGeometry> meshGeometry;
    std::shared_ptr<ConservationLaw> conservationLaw;
    std::shared_ptr<IntercellFluxScheme> intercellFluxScheme;
    std::shared_ptr<ConstrainedTransport> constrainedTransport;
    std::shared_ptr<BoundaryCondition> boundaryCondition;
    std::shared_ptr<RiemannSolver> riemannSolver;

    InitialDataFunction initialDataFunction;
    InitialDataFunction vectorPotentialFunction;
    InitialDataFunction boundaryValueFunction;
};

#endif
