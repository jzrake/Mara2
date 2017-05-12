#include <iostream>
#include "Configuration.hpp"
#include "CartesianMeshGeometry.hpp"
#include "sol.hpp"
#define lua (luaState->L)
#define PRINT_VEC3(v) "[" << v[0] << " " << v[1] << " " << v[2] << "]"




// ============================================================================
class Configuration::LuaState
{
public:
    sol::state L;

    InitialDataFunction makeIDF (sol::function& func)
    {
        return [=] (double x, double y, double z)
        {
            sol::table result = func (x, y, z);
            std::vector<double> P;

            for (int n = 0; n < result.size(); ++n)
            {
                P.push_back (result[n + 1]);
            }
            return P;
        };
    };
};




// ============================================================================
Configuration::Configuration()
{
    luaState.reset (new LuaState);
    lua.open_libraries (sol::lib::base, sol::lib::math);
}

Configuration::~Configuration()
{

}

SimulationSetup Configuration::fromLuaFile (std::string filename)
{
    std::cout << std::string(80, '-') << std::endl;
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    std::cout << std::string(80, '-') << std::endl;


    auto setup = SimulationSetup();
    lua.script_file (filename);

    sol::function initial_data = lua["initial_data"];
    setup.initialDataFunction = luaState->makeIDF (initial_data);


    // grid_geometry
    // ------------------------------------------------------------------------
    std::string grid_geometry = lua["grid_geometry"];
    std::cout << "grid geometry: " << grid_geometry << std::endl;

    if (grid_geometry == "cartesian")
    {
        auto geometry = new CartesianMeshGeometry;

        geometry->shape[0] = lua["resolution"][1].get_or (128);
        geometry->shape[1] = lua["resolution"][2].get_or (1);
        geometry->shape[2] = lua["resolution"][3].get_or (1);
        geometry->lower[0] = lua["domain_lower"][1].get_or (0.0);
        geometry->lower[1] = lua["domain_lower"][2].get_or (0.0);
        geometry->lower[2] = lua["domain_lower"][3].get_or (0.0);
        geometry->upper[0] = lua["domain_upper"][1].get_or (1.0);
        geometry->upper[1] = lua["domain_upper"][2].get_or (1.0);
        geometry->upper[2] = lua["domain_upper"][3].get_or (1.0);

        std::cout << "resolution: " << PRINT_VEC3(geometry->shape) << std::endl;
        std::cout << "domain_lower: " << PRINT_VEC3(geometry->lower) << std::endl;
        std::cout << "domain_upper: " << PRINT_VEC3(geometry->upper) << std::endl;

        setup.meshGeometry.reset (geometry);
    }
    else
    {
        throw std::runtime_error ("illegal option for grid_geometry");
    }


    // boundary condition
    // ------------------------------------------------------------------------
    std::string boundary_condition = lua["boundary_condition"];
    std::cout << "boundary_condition: " << boundary_condition << std::endl;

    if (boundary_condition == "periodic")
    {
        setup.boundaryCondition.reset (new PeriodicBoundaryCondition());
    }
    else
    {
        throw std::runtime_error ("illegal option for boundary_condition");
    }


    // conservation law
    // ------------------------------------------------------------------------
    std::string conservation_law = lua["conservation_law"][1];
    std::cout << "conservation_law: " << conservation_law << std::endl;

    if (conservation_law == "scalar_advection")
    {
        double wave_speed = lua["conservation_law"]["wave_speed"].get_or (1.0);
        std::cout << "wave_speed: " << wave_speed << std::endl;

        setup.conservationLaw.reset (new ScalarAdvection (wave_speed));
    }
    else
    {
        throw std::runtime_error ("illegal option for conservation_law");
    }


    // intercell flux scheme
    // ------------------------------------------------------------------------
    setup.intercellFluxScheme.reset (new ScalarUpwind);


    // Run configuration
    // ------------------------------------------------------------------------
    setup.outputDirectory = lua["output_directory"];
    setup.finalTime = lua["final_time"];
    setup.checkpointInterval = lua["checkpoint_interval"];
    setup.cflParameter = lua["cfl_parameter"];

    std::cout << "output_directory: " << setup.outputDirectory << std::endl;
    std::cout << "final_time: " << setup.finalTime << std::endl;
    std::cout << "checkpoint_interval: " << setup.checkpointInterval << std::endl;
    std::cout << "cfl_parameter: " << setup.cflParameter << std::endl;


    return setup;
}
