#include <iostream>
#include <iomanip>
#include "UserConfiguration.hpp"
#define GUARD_STRING std::string (48, '-')

UserConfiguration::UserConfiguration (int argc, const char* argv[])
{
    members["final_time"] = 0;
    members["checkpoint_interval"] = 1.0;
    members["vtk_output_interval"] = 1.0;
    members["cfl_parameter"] = 0.25;
    members["pressure_floor"] = 0.0;
    members["vtk_use_binary"] = true;
    members["disable_ct"] = true;
    members["output_directory"] = "data";
    members["run_name"] = "test";
    members["restart_file"] = "";
    members["num_cells_1"] = 128;
    members["num_cells_2"] = 1;
    members["num_cells_3"] = 1;
    members["domain_lower_1"] = -0.5;
    members["domain_lower_2"] = -0.5;
    members["domain_lower_3"] = -0.5;
    members["domain_upper_1"] = +0.5;
    members["domain_upper_2"] = +0.5;
    members["domain_upper_3"] = +0.5;
    members["runge_kutta_order"] = 0;
    members["flux_scheme"] = "method_of_lines_plm";
    members["boundary_condition"] = "periodic";
    members["mesh_geometry"] = "cartesian";
    members["riemann_solver"] = "hlle";
    members["conservation_law"] = "newtonian_hydro";
    members["constrained_transport"] = "uniform_cartesian";
    members["initial_data_function"] = "shocktube1";
    members["vector_potential_function"] = "";
    members["boundary_value_function"] = "";

    Variant::updateFromCommandLine (members, argc, argv);
}

const Variant::NamedValues& UserConfiguration::getMembers() const
{
    return members;
}

void UserConfiguration::describe() const
{
    std::cout << GUARD_STRING << std::endl;
    std::cout << "User configuration:" << std::endl;
    std::cout << GUARD_STRING << std::endl;
    std::cout << members;
    std::cout << GUARD_STRING << std::endl;
}
