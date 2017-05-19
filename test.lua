local setups = require 'setups'

 -- name of the current run
run_name = 'AdvectionTest'

-- directory where output should go
output_directory = 'data'

-- Time at which to terminate simulation
final_time = 4.0

-- How frequently to output full simulation snapshot
checkpoint_interval = 0.01

-- CFL parameter
cfl_parameter = 0.5

-- A callback to supply simple initial data
initial_data = setups['hydro_uniform']

-- Grid geometry
grid_geometry = 'cartesian'

-- Grid resolution (Must be a 3D array)
resolution = {256, 16, 1}

-- Domain lower bounds
domain_lower = { 0.0, 0.00,-0.125}

-- Domain upper bounds
domain_upper = { 1.0, 0.25, 0.125}

-- Conservation law: scalar_advection, newtonian_hydro, newtonian_mhd
conservation_law = {'newtonian_hydro'}

-- Riemann solver: upwind, hlle
riemann_solver = 'hlle'

-- Flux scheme: scalar_upwind, method_of_lines
flux_scheme = {'method_of_lines_plm', plm_theta=1.5}

-- RK order: must be 1, 2, or 3
runge_kutta_order = 2

-- Boundary condition name
boundary_condition = 'planar_pipe_flow'
