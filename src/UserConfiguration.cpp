#include <iostream>
#include <iomanip>
#include "UserConfiguration.hpp"


UserConfiguration::UserConfiguration()
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
    members["flux_scheme"] = "plm";
    members["conservation_law"] = "newtonian_hydro";
    members["initial_data_function"] = "shocktube1";
    members["vector_potential_function"] = "";
    members["boundary_value_function"] = "";
}

void UserConfiguration::describe() const
{
    std::cout << members;
}
