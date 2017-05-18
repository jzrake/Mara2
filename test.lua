local setups = require 'setups'

 -- name of the current run
run_name = 'AdvectionTest'

-- directory where output should go
output_directory = 'euler'

-- Time at which to terminate simulation
final_time = 0.1

-- How frequently to output full simulation snapshot
checkpoint_interval = 0.001

-- CFL parameter
cfl_parameter = 0.5

-- A callback to supply simple initial data
initial_data = setups['shocktube1']

-- Grid geometry
grid_geometry = 'cartesian'

-- Grid resolution (Must be a 3D array)
resolution = {300, 1, 1}

-- Domain lower bounds
domain_lower = {-1.0, 0.0, 0.0}

-- Domain upper bounds
domain_upper = { 1.0, 1.0, 1.0}

-- Fluid variables: scalar_advection, euler_equation
-- conservation_law = {'scalar_advection', wave_speed=1}
conservation_law = {'euler_equation'}

-- Riemann solver: upwind, hlle
riemann_solver = 'hlle'

-- Flux scheme: scalar_upwind, method_of_lines
flux_scheme = {'method_of_lines', plm_theta=2.0}

-- RK order: must be 1, 2, or 3
runge_kutta_order = 3

-- Boundary condition name
boundary_condition = 'periodic'
