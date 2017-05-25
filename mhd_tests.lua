local mara = require 'mara'

local fieldWaveNumber = 10.0
local fieldAmplitude = 1e-5


local function background(x, y, z)
	local v1 = 1.0
	local v2 = 0.0
	local v3 = 0.0
	return {1.0, v1, v2, v3, 1.0, 0.0, 0.0, 0.0}
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
	initial_data = background,
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
	setup.output_directory = test.output_directory
	setup.vector_potential = test.vector_potential
	mara.run(setup)
end


--
-- Do a single run in 2D advecting a field loop once around the domain
--
setup.resolution = {128, 128, 1}
setup.vector_potential = field_loop_z
setup.output_directory = 'field_loop_advect'
setup.final_time = 1.0
mara.run(setup)
