#include "SimulationSetup.hpp"
#include "CartesianMeshGeometry.hpp"
#include "BoundaryConditions.hpp"
#include "ConservationLaws.hpp"
#include "ConstrainedTransport.hpp"
#include "IntercellFluxSchemes.hpp"
#include "RiemannSolvers.hpp"







// ============================================================================
SimulationSetup2::SimulationSetup2 (Variant::NamedValues config)
{
	auto lookupBoundaryCondition = [] (std::string name) -> std::shared_ptr<BoundaryCondition>
	{
		if (name == "periodic")           return std::make_shared<PeriodicBoundaryCondition>();
		if (name == "reflecting")         return std::make_shared<ReflectingBoundaryCondition>();
		if (name == "outflow")            return std::make_shared<OutflowBoundaryCondition>();
		throw std::runtime_error ("no boundary condition named " + name);
	};
	auto lookupConservationLaw = [] (std::string name) -> std::shared_ptr<ConservationLaw>
	{
		if (name == "scalar_advection")   return std::make_shared<ScalarAdvection>();
		if (name == "newtonian_hydro")    return std::make_shared<NewtonianHydro>();
		if (name == "newtonian_mhd")      return std::make_shared<NewtonianMHD>();
		throw std::runtime_error ("no conservation law named " + name);
	};
	auto lookupMeshGeometryLaw = [] (std::string name) -> std::shared_ptr<MeshGeometry>
	{
		if (name == "cartesian")          return std::make_shared<CartesianMeshGeometry>();
		throw std::runtime_error ("no mesh geometry named " + name);
	};
	auto lookupFluxScheme = [] (std::string name) -> std::shared_ptr<IntercellFluxScheme>
	{
		if (name == "method_of_lines")        return std::make_shared<MethodOfLines>();
		if (name == "method_of_lines_plm")    return std::make_shared<MethodOfLinesPlm>();
		if (name == "method_of_lines_weno")   return std::make_shared<MethodOfLinesWeno>();
		throw std::runtime_error ("no mesh geometry named " + name);
	};
	auto lookupRiemannSolver = [] (std::string name) -> std::shared_ptr<RiemannSolver>
	{
		if (name == "upwind")             return std::make_shared<UpwindRiemannSolver>();
		if (name == "hlle")               return std::make_shared<HlleRiemannSolver>();
		throw std::runtime_error ("no riemann solver named " + name);
	};
	auto lookupConstrainedTransport = [] (std::string name) -> std::shared_ptr<ConstrainedTransport>
	{
		if (name == "uniform_cartesian")  return std::make_shared<UniformCartesianCT>();
		throw std::runtime_error ("no constrained transport named " + name);
	};
	auto lookupInitialDataFunction = [] (std::string name) -> InitialDataFunction
	{
		return nullptr;
	};

	boundaryCondition       = lookupBoundaryCondition    (config["boundary_condition"]);
	conservationLaw         = lookupConservationLaw      (config["conservation_law"]);
	meshGeometry            = lookupMeshGeometryLaw      (config["mesh_geometry"]);
	riemannSolver           = lookupRiemannSolver        (config["riemann_solver"]);
	intercellFluxScheme     = lookupFluxScheme           (config["flux_scheme"]);
	constrainedTransport    = lookupConstrainedTransport (config["constrained_transport"]);
	initialDataFunction     = lookupInitialDataFunction  (config["initial_data_function"]);
	vectorPotentialFunction = lookupInitialDataFunction  (config["vector_potential_function"]);
	boundaryValueFunction   = lookupInitialDataFunction  (config["boundary_value_function"]);
}
