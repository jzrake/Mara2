#ifndef UserConfiguration_hpp
#define UserConfiguration_hpp

#include "Variant.hpp"




class UserConfiguration
{
public:
    UserConfiguration();
    void describe() const;

    bool vtk_use_binary;
    bool disable_ct;
    double final_time;
    double checkpoint_interval;
    double vtk_output_interval;
    double cfl_parameter;
    double plm_theta;
    double scalar_wavespeed;
    double adiabatic_gamma;
    double pressure_floor;
    int runge_kutta_order;
    int num_cells_1;
    int num_cells_2;
    int num_cells_3;
    double domain_lower_1;
    double domain_lower_2;
    double domain_lower_3;
    double domain_upper_1;
    double domain_upper_2;
    double domain_upper_3;
    std::string output_directory;
    std::string run_name;
    std::string restart_file;   
    std::string boundary_condition;
    std::string riemann_solver;
    std::string flux_scheme;
    std::string initial_data_function;
    std::string vector_potential_function;
    std::string boundary_value_function;

private:
    Variant::NamedValues members;
};

#endif
