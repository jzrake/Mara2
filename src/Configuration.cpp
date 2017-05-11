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


    // output options
    // ------------------------------------------------------------------------
    std::string outdir = lua["outdir"];
    std::cout << "outdir: " << outdir << std::endl;


    // final time, etc
    // ------------------------------------------------------------------------
    double final_time = lua["final_time"];
    double checkpoint_interval = lua["checkpoint_interval"];
    double cfl_parameter = lua["cfl_parameter"];
    std::cout << "final_time: " << final_time << std::endl;
    std::cout << "checkpoint_interval: " << checkpoint_interval << std::endl;
    std::cout << "cfl_parameter: " << cfl_parameter << std::endl;

    setup.finalTime = final_time;
    setup.checkpointInterval = checkpoint_interval;
    setup.cflParameter = cfl_parameter;


    setup.conservationLaw.reset (new ScalarAdvection);
    setup.intercellFluxScheme.reset (new ScalarUpwind);


    return setup;
}
