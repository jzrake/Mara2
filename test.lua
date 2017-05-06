
 -- name of the current run
run_name = 'AdvectionTest'

-- directory where output should go
outdir = './data'

-- Grid resolution (Must be a 3D array)
resolution = {128, 1, 1}

N1 = {123}

-- Fluid variables
conservation_law = 'scalar_advection'

-- Name of scheme
scheme = "method_of_lines"

-- CFL parameter
cfl_parameter = 0.5

initial_data = function() return {99, 33, 22} end
