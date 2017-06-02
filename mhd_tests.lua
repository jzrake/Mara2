local mara = require 'mara'

local fieldWaveNumber = 10.0
local fieldAmplitude = 1e-5


local function background(x, y, z)
	return {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0}
end

local function background_x(x, y, z)
	return {1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0}
end

local function background_y(x, y, z)
	return {1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0}
end

local function background_z(x, y, z)
	return {1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0}
end

local function field_loop_x(x, y, z)
	local R = (y^2 + z^2)^0.5
	local A = fieldAmplitude * math.exp (-(fieldWaveNumber * R)^3)
	return {A, 0.0, 0.0}
end

local function field_loop_y(x, y, z)
	local R = (z^2 + x^2)^0.5
	local A = fieldAmplitude * math.exp (-(fieldWaveNumber * R)^3)
	return {0.0, A, 0.0}
end

local function field_loop_z(x, y, z)
	local R = (x^2 + y^2)^0.5
	local A = fieldAmplitude * math.exp (-(fieldWaveNumber * R)^3)
	return {0.0, 0.0, A}
end

local function abc_equilibrium(x, y, z)
	local k = 4.0 * math.pi
	local A = 0.1 / k
	local B = 0.1 / k
	local C = 0.1 / k
	local a1 = C * math.cos(k * z) - B * math.sin(k * y)
	local a2 = A * math.cos(k * x) - C * math.sin(k * z)
	local a3 = B * math.cos(k * y) - A * math.sin(k * x)
	return {a1, a2, a3}
end

local function new_setup(args)
	local setup = {
		run_name = args['run_name'],
		output_directory = 'data/' .. args['run_name'],
		final_time = 1.0,
		checkpoint_interval = 0.1,
		vtk_output_interval = 0.025,
		vtk_use_binary = true,
		cfl_parameter = 0.3,
		grid_geometry = 'cartesian',
		resolution = {32, 32, 32},
		domain_lower = {-0.5, -0.5, -0.5},
		domain_upper = { 0.5,  0.5,  0.5},
		conservation_law = {'newtonian_mhd'},
		riemann_solver = 'hlle',
		flux_scheme = {'method_of_lines_plm', plm_theta=1.5},
		runge_kutta_order = 2,
		boundary_condition = 'periodic',
	}
	for k, v in pairs(args) do
		setup[k] = v
	end
	return setup
end

local tests = { }

--
-- Generate initial data for field loops around each axis. No evolution is
-- done here, we just inspect the output files to ensure the magnetic field
-- is divergenceless.
--
tests[1] = new_setup {run_name='field_loop_x', final_time=0.0, initial_data=background, vector_potential=field_loop_x}
tests[2] = new_setup {run_name='field_loop_y', final_time=0.0, initial_data=background, vector_potential=field_loop_y}
tests[3] = new_setup {run_name='field_loop_z', final_time=0.0, initial_data=background, vector_potential=field_loop_z}

--
-- Run a 2D advecting a field loop once around the domain, with
-- the motion along along one axis in the plane. There is one run
-- for each of the three axes, so we test all of the CT algorithm.
--
tests[4] = new_setup {run_name='field_loop_advect_x', resolution={1, 64, 64}, initial_data=background_x, vector_potential=field_loop_x}
tests[5] = new_setup {run_name='field_loop_advect_y', resolution={64, 1, 64}, initial_data=background_y, vector_potential=field_loop_y}
tests[6] = new_setup {run_name='field_loop_advect_z', resolution={64, 64, 1}, initial_data=background_z, vector_potential=field_loop_z}

--
-- These are 3D magnetic equilibrium tests.
--
tests[7] = new_setup {run_name='abc_equilibrium', initial_data=background, vector_potential=abc_equilibrium}

for _, setup in pairs (tests) do
	mara.run(setup)
end
