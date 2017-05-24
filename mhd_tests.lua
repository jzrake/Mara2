--local mara = require 'mara'
local setups = require 'setups'

local fieldLoopTest = {
	run_name = 'FieldLoopTest',
	output_directory = 'data/field_loop_test',
	final_time = 0.0,
	checkpoint_interval = 0.1,
	vtk_output_interval = 0.01,
	cfl_parameter = 0.3,
	initial_data = setups['mhd_field_loop'],
	vector_potential = setups['mhd_field_loop_A'],
	grid_geometry = 'cartesian',
	resolution = {16, 16, 1},
	domain_lower = {-0.5, -0.5, -0.5},
	domain_upper = { 0.5,  0.5,  0.5},
	conservation_law = {'newtonian_mhd'},
	riemann_solver = 'hlle',
	flux_scheme = {'method_of_lines_plm', plm_theta=1.5},
	runge_kutta_order = 2,
	boundary_condition = 'periodic',
}

--mara.run (fieldLoopTest)


for k, v in pairs (fieldLoopTest) do
	_G[k] = v
end

