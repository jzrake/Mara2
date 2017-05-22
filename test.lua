local setups = require 'setups'

-- name of the current run
run_name = 'AdvectionTest'

-- directory where output should go
output_directory = 'data'

-- Time at which to terminate simulation
final_time = 16.0

-- How frequently to output full simulation snapshot (set to 0 to disable
-- checkpoints)
checkpoint_interval = 1.0

-- How frequently to write VTK data (set to 0 to disable VTK output)
vtk_output_interval = 0.05

-- CFL parameter
cfl_parameter = 0.3

-- A callback to supply simple initial data
initial_data = setups['mhd_uniform_vertical']

-- Grid geometry
grid_geometry = 'cartesian'

-- Grid resolution (Must be a 3D array)
resolution = {32, 32, 32}

-- Domain lower bounds
domain_lower = {-0.5, -0.5, -0.5}

-- Domain upper bounds
domain_upper = { 0.5,  0.5,  0.5}

-- Conservation law: scalar_advection, newtonian_hydro, newtonian_mhd
conservation_law = {'newtonian_mhd'}

-- Riemann solver: upwind, hlle
riemann_solver = 'hlle'

-- Flux scheme: scalar_upwind, method_of_lines
flux_scheme = {'method_of_lines_plm', plm_theta=1.5}

-- RK order: must be 1, 2, or 3
runge_kutta_order = 2

-- Boundary condition name
boundary_condition = 'driven_mhd'
