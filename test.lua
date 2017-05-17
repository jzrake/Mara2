
 -- name of the current run
run_name = 'AdvectionTest'

-- directory where output should go
output_directory = './'

-- Time at which to terminate simulation
final_time = 0.1

-- How frequently to output full simulation snapshot
checkpoint_interval = 1.0

-- CFL parameter
cfl_parameter = 0.5

-- A callback to supply simple initial data
-- initial_data = function(x, y, z) return {math.exp(-x^2 / 0.025), 2 - x*x} end

initial_data = function(x, y, z)
	local d = 1.0 + math.exp(-x^2 / 0.025)
	local u = 0.0
	local v = 0.0
	local w = 0.0
	local p = 1.0
	return {d, u, v, w, p}
end

-- Grid geometry
grid_geometry = 'cartesian'

-- Grid resolution (Must be a 3D array)
resolution = {100, 1, 1}

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
