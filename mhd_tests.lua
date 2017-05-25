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

local setup = {
	run_name = 'FieldLoopTest',
	final_time = 0.0,
	checkpoint_interval = 0.0,
	vtk_output_interval = 0.025,
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




--
-- Generate initial data for field loops around each axis. No evolution is
-- done here, we just inspect the output files to ensure the magnetic field
-- is divergenceless.
--
local tests = {
	{vector_potential = field_loop_x, output_directory = 'field_loop_x'},
	{vector_potential = field_loop_y, output_directory = 'field_loop_y'},
	{vector_potential = field_loop_z, output_directory = 'field_loop_z'},
}

for _, test in ipairs(tests) do
	setup.initial_data = background_x
	setup.vector_potential = test.vector_potential
	setup.output_directory = test.output_directory
	mara.run(setup)
end



--
-- Run a 2D advecting a field loop once around the domain, with
-- the motion along along one axis in the plane. There is one run
-- for each of the three axes, so we test all of the CT algorithm.
--
setup.final_time = 1.0

setup.resolution = {1, 64, 64}
setup.initial_data = background_y
setup.vector_potential = field_loop_x
setup.output_directory = 'field_loop_advect_x'
mara.run(setup)

setup.resolution = {64, 1, 64}
setup.initial_data = background_z
setup.vector_potential = field_loop_y
setup.output_directory = 'field_loop_advect_y'
mara.run(setup)

setup.resolution = {64, 64, 1}
setup.initial_data = background_x
setup.vector_potential = field_loop_z
setup.output_directory = 'field_loop_advect_z'
mara.run(setup)




--
-- These are 3D magnetic equilibrium tests.
--
local function abc_equilibrium_A(x, y, z)
	local k = 4.0 * math.pi
	local A = 0.1 / k
	local B = 0.1 / k
	local C = 0.1 / k
	local a1 = C * math.cos(k * z) - B * math.sin(k * y)
	local a2 = A * math.cos(k * x) - C * math.sin(k * z)
	local a3 = B * math.cos(k * y) - A * math.sin(k * x)
	return {a1, a2, a3}
end
setup.final_time = 1.0
setup.resolution = {32, 32, 1}
setup.initial_data = background
setup.vector_potential = abc_equilibrium_A
setup.output_directory = 'abc_equilibrium'
mara.run(setup)
