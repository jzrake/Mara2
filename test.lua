

 -- name of the current run
run_name = 'AdvectionTest'

-- directory where output should go
output_directory = './'

-- Time at which to terminate simulation
final_time = 0.5

-- How frequently to output full simulation snapshot
checkpoint_interval = 1.0

-- CFL parameter
cfl_parameter = 0.5

-- A callback to supply simple initial data
initial_data = function(x, y, z) return {math.exp(-x^2 / 0.025)} end

-- Grid geometry
grid_geometry = 'cartesian'

-- Grid resolution (Must be a 3D array)
resolution = {128, 1, 1}

-- Domain lower bounds
domain_lower = {-1.0, 0.0, 0.0}

-- Domain upper bounds
domain_upper = { 1.0, 1.0, 1.0}

-- Fluid variables
conservation_law = {'scalar_advection', wave_speed=1}

-- Boundary condition name
boundary_condition = 'periodic'

-- Name of scheme
scheme = "method_of_lines"
