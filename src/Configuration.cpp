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

    InitialDataFunction makeIDF (sol::function func);
    Cow::Shape makeShape (sol::table table);
    MeshGeometry::Coordinate makeCoordinate (sol::table table);
    SimulationSetup fromLuaTable (sol::table);
};




// ============================================================================
InitialDataFunction Configuration::LuaState::makeIDF (sol::function func)
{
    if (func == sol::nil)
    {
        return nullptr;
    }
    return [=] (double x, double y, double z)
    {
        sol::table result = func (x, y, z);
        std::vector<double> P;

        for (unsigned int n = 0; n < result.size(); ++n)
        {
            P.push_back (result[n + 1]);
        }
        return P;
    };
};

Cow::Shape Configuration::LuaState::makeShape (sol::table table)
{
    if (table == sol::nil)
    {
        throw std::runtime_error ("[Configuration]: expected 'shape' like {x, y, z}, got nil");
    }
    auto shape = Cow::Shape();
    shape[0] = table[1];
    shape[1] = table[2];
    shape[2] = table[3];
    shape[3] = 1;
    shape[4] = 1;
    return shape;
};

MeshGeometry::Coordinate Configuration::LuaState::makeCoordinate (sol::table table)
{
    if (table == sol::nil)
    {
        throw std::runtime_error ("[Configuration]: expected 'coordiante' like {x, y, z}, got nil");
    }
    auto X = MeshGeometry::Coordinate();
    X[0] = table[1];
    X[1] = table[2];
    X[2] = table[3];
    return X;
}


SimulationSetup Configuration::LuaState::fromLuaTable (sol::table cfg)
{
    auto lookupBoundaryCondition = [] (std::string name) -> std::shared_ptr<BoundaryCondition>
    {
        if (name == "periodic")               return std::make_shared<PeriodicBoundaryCondition>();
        if (name == "reflecting")             return std::make_shared<ReflectingBoundaryCondition>();
        if (name == "outflow")                return std::make_shared<OutflowBoundaryCondition>();
        if (name == "driven_mhd")             return std::make_shared<DrivenMHDBoundary>();
        if (name.empty()) throw std::runtime_error ("no option given for boundary_condition");
        throw std::runtime_error ("no boundary condition named " + name);
    };
    auto lookupConservationLaw = [] (std::string name) -> std::shared_ptr<ConservationLaw>
    {
        if (name == "scalar_advection")       return std::make_shared<ScalarAdvection>();
        if (name == "newtonian_hydro")        return std::make_shared<NewtonianHydro>();
        if (name == "newtonian_mhd")          return std::make_shared<NewtonianMHD>();
        if (name.empty()) throw std::runtime_error ("no option given for conservation_law");
        throw std::runtime_error ("no conservation law named " + name);
    };
    auto lookupMeshGeometryLaw = [] (std::string name) -> std::shared_ptr<MeshGeometry>
    {
        if (name == "cartesian")              return std::make_shared<CartesianMeshGeometry>();
        if (name.empty()) throw std::runtime_error ("no option given for mesh_geometry");
        throw std::runtime_error ("no mesh geometry named " + name);
    };
    auto lookupRiemannSolver = [] (std::string name) -> std::shared_ptr<RiemannSolver>
    {
        if (name == "upwind")                 return std::make_shared<UpwindRiemannSolver>();
        if (name == "hlle")                   return std::make_shared<HlleRiemannSolver>();
        if (name.empty()) throw std::runtime_error ("no option given for riemann_solver");
        throw std::runtime_error ("no riemann solver named " + name);
    };
    auto lookupFluxScheme = [] (std::string name) -> std::shared_ptr<IntercellFluxScheme>
    {
        if (name == "method_of_lines")        return std::make_shared<MethodOfLines>();
        if (name == "method_of_lines_plm")    return std::make_shared<MethodOfLinesPlm>();
        if (name == "method_of_lines_weno")   return std::make_shared<MethodOfLinesWeno>();
        if (name.empty()) throw std::runtime_error ("no option given for flux_scheme");
        throw std::runtime_error ("no flux scheme named " + name);
    };
    auto lookupConstrainedTransport = [] (std::string name) -> std::shared_ptr<ConstrainedTransport>
    {
        if (name == "uniform_cartesian")      return std::make_shared<UniformCartesianCT>();
        if (name.empty()) throw std::runtime_error ("no option given for constrained_transport");
        throw std::runtime_error ("no constrained transport named " + name);
    };

    auto setup = SimulationSetup();
    setup.boundaryCondition       = lookupBoundaryCondition    (cfg["boundary_condition"]);
    setup.conservationLaw         = lookupConservationLaw      (cfg["conservation_law"]);
    setup.meshGeometry            = lookupMeshGeometryLaw      (cfg["mesh_geometry"]);
    setup.riemannSolver           = lookupRiemannSolver        (cfg["riemann_solver"]);
    setup.intercellFluxScheme     = lookupFluxScheme           (cfg["flux_scheme"]);
    setup.constrainedTransport    = lookupConstrainedTransport (cfg["constrained_transport"]);
    setup.initialDataFunction     = makeIDF                    (cfg["initial_data"]);
    setup.vectorPotentialFunction = makeIDF                    (cfg["vector_potential"]);
    setup.boundaryValueFunction   = makeIDF                    (cfg["boundary_value"]);
    setup.finalTime               =                             cfg["final_time"];
    setup.outputDirectory         =                             cfg["output_directory"];
    setup.checkpointInterval      =                             cfg["checkpoint_interval"];
    setup.vtkOutputInterval       =                             cfg["vtk_output_interval"];
    setup.timeSeriesInterval      =                             cfg["time_series_interval"];
    setup.cflParameter            =                             cfg["cfl_parameter"];
    setup.rungeKuttaOrder         =                             cfg["runge_kutta_order"].get_or (2);
    setup.disableCT               =                             cfg["disable_ct"].get_or (false);
    setup.runName                 =                             cfg["run_name"];
    setup.vtkUseBinary            =                             cfg["vtk_use_binary"].get_or (true);

    setup.conservationLaw->setAdvectionSpeed (cfg["wave_speed"].get_or (1.0), 0.0, 0.0);
    setup.conservationLaw->setPressureFloor (cfg["pressure_floor"].get_or (-1.0));
    setup.meshGeometry->setCellsShape (makeShape (cfg["resolution"]));
    setup.meshGeometry->setLowerUpper (
        makeCoordinate (cfg["domain_lower"]),
        makeCoordinate (cfg["domain_upper"]));

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
