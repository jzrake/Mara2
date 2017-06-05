#include <iostream>
#include <fstream>
#include <sstream>
#include "Configuration.hpp"
#include "CartesianMeshGeometry.hpp"
#include "BoundaryConditions.hpp"
#include "ConservationLaws.hpp"
#include "ConstrainedTransport.hpp"
#include "IntercellFluxSchemes.hpp"
#include "RiemannSolvers.hpp"

// Cow
#include "HDF5.hpp"

// Lua
#define lua (luaState->L)
// #define SOL_CHECK_ARGUMENTS // <-- causes stack error?
#define SOL_SAFE_USERTYPES
#define SOL_SAFE_FUNCTION
#include "sol.hpp"




class Configuration::LuaState
{
public:
    sol::state L;

    InitialDataFunction makeIDF (sol::function& func);
    SimulationSetup fromLuaTable (sol::table);
};




// ============================================================================
InitialDataFunction Configuration::LuaState::makeIDF (sol::function& func)
{
    if (func == sol::nil)
    {
        return nullptr;
    }
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

SimulationSetup Configuration::LuaState::fromLuaTable (sol::table cfg)
{
    auto setup = SimulationSetup();


    // initial data function
    // ------------------------------------------------------------------------
    sol::function initial_data = cfg["initial_data"];
    sol::function vector_potential = cfg["vector_potential"];
    setup.initialDataFunction = makeIDF (initial_data);
    setup.vectorPotentialFunction = makeIDF (vector_potential);


    // grid geometry
    // ------------------------------------------------------------------------
    std::string grid_geometry = cfg["grid_geometry"];

    if (grid_geometry == "cartesian")
    {
        auto geometry = new CartesianMeshGeometry;

        geometry->shape[0] = cfg["resolution"][1].get_or (128);
        geometry->shape[1] = cfg["resolution"][2].get_or (1);
        geometry->shape[2] = cfg["resolution"][3].get_or (1);
        geometry->lower[0] = cfg["domain_lower"][1].get_or (0.0);
        geometry->lower[1] = cfg["domain_lower"][2].get_or (0.0);
        geometry->lower[2] = cfg["domain_lower"][3].get_or (0.0);
        geometry->upper[0] = cfg["domain_upper"][1].get_or (1.0);
        geometry->upper[1] = cfg["domain_upper"][2].get_or (1.0);
        geometry->upper[2] = cfg["domain_upper"][3].get_or (1.0);

        setup.meshGeometry.reset (geometry);
    }
    else
    {
        throw std::runtime_error ("unrecognized option for grid_geometry");
    }


    // boundary condition
    // ------------------------------------------------------------------------
    std::string boundary_condition = cfg["boundary_condition"];

    if (boundary_condition == "periodic")
    {
        setup.boundaryCondition.reset (new PeriodicBoundaryCondition);
    }
    else if (boundary_condition == "outflow")
    {
        setup.boundaryCondition.reset (new OutflowBoundaryCondition);
    }
    else if (boundary_condition == "reflecting")
    {
        setup.boundaryCondition.reset (new ReflectingBoundaryCondition);
    }
    else if (boundary_condition == "driven_mhd")
    {
        auto drvBoundary = new DrivenMHDBoundary;
        sol::function F = cfg["boundary_velocity_function"];
        drvBoundary->setBoundaryValueFunction (makeIDF (F));
        setup.boundaryCondition.reset (drvBoundary);
    }
    else
    {
        throw std::runtime_error ("unrecognized option for boundary_condition");
    }


    // conservation law
    // ------------------------------------------------------------------------
    std::string conservation_law = cfg["conservation_law"][1];

    if (conservation_law == "scalar_advection")
    {
        double wave_speed = cfg["conservation_law"]["wave_speed"].get_or (1.0);
        setup.conservationLaw.reset (new ScalarAdvection);
        setup.conservationLaw->setAdvectionSpeed (wave_speed, 0.0, 0.0);
    }
    else if (conservation_law == "newtonian_hydro")
    {
        setup.conservationLaw.reset (new NewtonianHydro);
    }
    else if (conservation_law == "newtonian_mhd")
    {
        setup.conservationLaw.reset (new NewtonianMHD);
    }
    else
    {
        throw std::runtime_error ("unrecognized option for conservation_law");
    }

    if (cfg["pressure_floor"])
    {
        setup.conservationLaw->setPressureFloor (cfg["pressure_floor"]);
    }


    // riemann solver
    // ------------------------------------------------------------------------
    std::string riemann_solver = cfg["riemann_solver"];

    if (riemann_solver == "upwind")
    {
        setup.riemannSolver.reset (new UpwindRiemannSolver);
    }
    else if (riemann_solver == "hlle")
    {
        setup.riemannSolver.reset (new HlleRiemannSolver);
    }
    else
    {
        throw std::runtime_error ("unrecognized option for riemann_solver");
    }


    // intercell flux scheme
    // ------------------------------------------------------------------------
    std::string flux_scheme = cfg["flux_scheme"][1];

    if (flux_scheme == "scalar_upwind")
    {
        setup.intercellFluxScheme.reset (new ScalarUpwind);
    }
    else if (flux_scheme == "method_of_lines")
    {
        auto scheme = new MethodOfLines;
        scheme->setRiemannSolver (setup.riemannSolver);
        setup.intercellFluxScheme.reset (scheme);
    }
    else if (flux_scheme == "method_of_lines_plm")
    {
        double plm_theta = cfg["flux_scheme"]["plm_theta"].get_or (1.5);
        auto scheme = new MethodOfLinesPlm (plm_theta);
        scheme->setRiemannSolver (setup.riemannSolver);
        setup.intercellFluxScheme.reset (scheme);
    }
    else if (flux_scheme == "method_of_lines_weno")
    {
        setup.intercellFluxScheme.reset (new MethodOfLinesWeno);
    }
    else
    {
        throw std::runtime_error ("unrecognized option for flux_scheme");
    }


    // constrained transport algorithm
    // ------------------------------------------------------------------------
    setup.constrainedTransport.reset (new UniformCartesianCT);


    // Run configuration
    // ------------------------------------------------------------------------
    setup.finalTime = cfg["final_time"];
    setup.outputDirectory = cfg["output_directory"];
    setup.checkpointInterval = cfg["checkpoint_interval"];
    setup.vtkOutputInterval = cfg["vtk_output_interval"];
    setup.cflParameter = cfg["cfl_parameter"];
    setup.rungeKuttaOrder = cfg["runge_kutta_order"];
    setup.disableCT = cfg["disable_ct"];
    setup.runName = cfg["run_name"];
    setup.vtkUseBinary = cfg["vtk_use_binary"].get_or (true);

    if (! (
        setup.rungeKuttaOrder == 1
     || setup.rungeKuttaOrder == 2
     || setup.rungeKuttaOrder == 3))
    {
        throw std::runtime_error ("runge_kutta_order must be 1, 2, or 3");
    }

    return setup;
}




// ============================================================================
Configuration::Configuration (int argc, const char* argv[])
{
    luaState.reset (new LuaState);
    lua.open_libraries (sol::lib::base, sol::lib::math, sol::lib::package);

    auto luaTableStream = std::ostringstream();
    bool passedDoubleDash = false;

    for (int n = 0; n < argc; ++n)
    {
        if (argv[n] == std::string ("--"))
        {
            passedDoubleDash = true;
            continue;
        }
        if (passedDoubleDash)
        {
            luaTableStream << argv[n] << "; ";
        }
    }
    commandLineLuaString = luaTableStream.str();
}

Configuration::~Configuration()
{

}

SimulationSetup Configuration::fromCheckpoint (std::string filename)
{
    auto file = Cow::H5::File (filename, "r");
    auto script = file.readString ("script");
    lua.script (script);
    lua.script (commandLineLuaString);
    auto setup = luaState->fromLuaTable (lua.globals());
    setup.luaScript = script;
    setup.restartFile = filename;
    setup.initialDataFunction = nullptr;
    setup.vectorPotentialFunction = nullptr;
    return setup;
}

SimulationSetup Configuration::fromLuaFile (std::string filename)
{
    lua.script_file (filename);
    lua.script (commandLineLuaString);
    auto setup = luaState->fromLuaTable (lua.globals());
    auto stream = std::ifstream (filename);
    auto c0 = std::istreambuf_iterator<char>(stream);
    auto c1 = std::istreambuf_iterator<char>();
    setup.luaScript = std::string (c0, c1);
    return setup;
}

int Configuration::launchFromScript (MaraSession& session, std::string filename)
{
    if (! commandLineLuaString.empty())
    {
        std::cout << "[Configuration] Warning: command line string ignored in run mode" << std::endl;
        std::cout << "[Configuration] " << commandLineLuaString << std::endl;
    }

    sol::table mara = lua.require_script ("mara", "return {}", false);
    mara["run"] = [&] (sol::table T)
    {
        auto setup = luaState->fromLuaTable (T);
        return session.launch (setup);
    };

    lua.script_file (filename);
    return 0;
}
