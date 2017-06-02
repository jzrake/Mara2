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
        drvBoundary->setVelocityFunction (makeIDF (F));
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
        setup.conservationLaw.reset (new ScalarAdvection (wave_speed));
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





// ============================================================================
#define NO_EXPERIMENT_CFG
#ifndef NO_EXPERIMENT_CFG

#include <iomanip>

#define SCL(nm) cout << left << setw(W) << setfill('.') << #nm << " " << nm << endl
#define VCI(nm) cout << left << setw(W) << setfill('.') << #nm << " " << nm[0] << " " << nm[1] << " " << nm[2] << " " << endl
#define VCD(nm) cout << left << setw(W) << setfill('.') << #nm << " " << showpos << nm[0] << " " << nm[1] << " " << nm[2] << " " << endl
#define FUN(nm) cout << left << setw(W) << setfill('.') << #nm << " " << (nm ? "f(...)" : "nil") << std::endl;
#define MEM(nm) #nm, &UserConfiguration::nm

namespace sol
{
    template <> struct lua_size<std::array<int, 3>> : std::integral_constant<int, 1> {};
    template <> struct lua_type_of<std::array<int, 3>> : std::integral_constant<sol::type, sol::type::table> {};

    namespace stack
    {
        template <>
        struct checker<std::array<int, 3>>
        {
            template <typename Handler>
            static bool check(lua_State* L, int index, Handler&& handler, record& tracking)
            {
                tracking.use(1);
                return stack::check<sol::table> (L, index, handler);
            }
        };

        template <>
        struct getter<std::array<int, 3>>
        {
            static std::array<int, 3> get(lua_State* L, int index, record& tracking)
            {
                sol::table tab = stack::get<sol::table> (L, index);
                tracking.use(1);
                return {{tab[1], tab[2], tab[3]}};
            }
        };

        template <>
        struct pusher<std::array<int, 3>>
        {
            static int push(lua_State* L, const std::array<int, 3>& arr)
            {
                return stack::push (L, sol::as_table (arr));
            }
        };

        template <>
        struct checker<std::array<double, 3>>
        {
            template <typename Handler>
            static bool check(lua_State* L, int index, Handler&& handler, record& tracking)
            {
                tracking.use(1);
                return stack::check<sol::table> (L, index, handler);
            }
        };

        template <>
        struct getter<std::array<double, 3>>
        {
            static std::array<double, 3> get(lua_State* L, int index, record& tracking)
            {
                sol::table tab = stack::get<sol::table> (L, index);
                tracking.use(1);
                return {{tab[1], tab[2], tab[3]}};
            }
        };

        template <>
        struct pusher<std::array<double, 3>>
        {
            static int push(lua_State* L, const std::array<double, 3>& arr)
            {
                return stack::push (L, sol::as_table (arr));
            }
        };
    }
}

class UserConfiguration
{
public:
    double final_time;
    double checkpoint_interval;
    double vtk_output_interval;
    double cfl_parameter;
    double pressure_floor;
    int runge_kutta_order;
    bool vtk_use_binary;
    bool disable_ct;
    std::string output_directory;
    std::string run_name;
    std::string lua_script;
    std::string restart_file;
    std::array<int, 3> domain_shape;
    std::array<double, 3> domain_lower;
    std::array<double, 3> domain_upper;
    InitialDataFunction initial_data_function;
    InitialDataFunction vector_potential_function;
    InitialDataFunction boundary_value_function;

    UserConfiguration()
    {
        final_time = 0;
        checkpoint_interval = 1.0;
        vtk_output_interval = 1.0;
        cfl_parameter = 0.25;
        pressure_floor = 0.0;
        runge_kutta_order = 0;
        vtk_use_binary = true;
        disable_ct = false;
        output_directory = "data";
        run_name = "test";
        lua_script = "";
        restart_file = "";
        domain_shape = {{ 128, 1, 1 }};
        domain_lower = {{-0.5,-0.5,-0.5 }};
        domain_upper = {{ 0.5, 0.5, 0.5 }};
        initial_data_function = nullptr;
        vector_potential_function = nullptr;
        boundary_value_function = nullptr;
    }

    void describe() const
    {
        using std::cout;
        using std::endl;
        using std::left;
        using std::setw;
        using std::setfill;
        using std::showpos;
        const int W = 32;

        SCL(final_time);
        SCL(checkpoint_interval);
        SCL(vtk_output_interval);
        SCL(cfl_parameter);
        SCL(pressure_floor);
        SCL(runge_kutta_order);
        SCL(vtk_use_binary);
        SCL(disable_ct);
        SCL(output_directory);
        SCL(run_name);
        SCL(lua_script);
        SCL(restart_file);
        VCI(domain_shape);
        VCD(domain_lower);
        VCD(domain_upper);
        FUN(initial_data_function);
        FUN(vector_potential_function);
        FUN(boundary_value_function);
    }
};

int Configuration::experiment (MaraSession& session, std::string filename)
{
    sol::table mara = lua.require_script ("mara", "return {}", false);
    mara.new_simple_usertype<UserConfiguration> ("Configuration",
        MEM(final_time),
        MEM(checkpoint_interval),
        MEM(vtk_output_interval),
        MEM(cfl_parameter),
        MEM(pressure_floor),
        MEM(runge_kutta_order),
        MEM(vtk_use_binary),
        MEM(disable_ct),
        MEM(output_directory),
        MEM(run_name),
        MEM(lua_script),
        MEM(restart_file),
        MEM(domain_shape),
        MEM(domain_lower),
        MEM(domain_upper),
        MEM(initial_data_function),
        MEM(vector_potential_function),
        MEM(boundary_value_function),
        MEM(describe));

    mara["run"] = [&] (sol::table T)
    {
        auto setup = luaState->fromLuaTable (T);
        return session.launch (setup);
    };

    lua.script_file (filename);
    return 0;
}

#else
int Configuration::experiment (MaraSession& session, std::string filename)
{
    std::cout << "experimental config is disabled";
    return 0;
}
#endif
