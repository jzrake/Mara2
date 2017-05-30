local setups = require 'setups'

-- name of the current run
run_name = 'ShocktubeTest'

-- directory where output should go
output_directory = 'test'

-- Time at which to terminate simulation
final_time = 0.0

-- How frequently to output full simulation snapshot (set to 0 to disable
-- checkpoints)
checkpoint_interval = 0.0

-- How frequently to write VTK data (set to 0 to disable VTK output)
vtk_output_interval = 0.0

-- CFL parameter
cfl_parameter = 0.6

-- A callback to supply simple initial data
initial_data = setups['shocktube1']

-- Grid geometry
grid_geometry = 'cartesian'

-- Grid resolution (Must be a 3D array)
resolution = {256, 1, 1}

-- Domain lower bounds
domain_lower = {-0.5, -0.5, -0.5}

-- Domain upper bounds
domain_upper = { 0.5,  0.5,  0.5}

-- Conservation law: scalar_advection, newtonian_hydro, newtonian_mhd
conservation_law = {'newtonian_hydro'}

-- Riemann solver: upwind, hlle
riemann_solver = 'hlle'

-- Flux scheme: scalar_upwind, method_of_lines
flux_scheme = {'method_of_lines_plm', plm_theta=2.0}

-- RK order: must be 1, 2, or 3
runge_kutta_order = 2

-- Boundary condition name
boundary_condition = 'periodic'
