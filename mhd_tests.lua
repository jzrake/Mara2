local mara = require 'mara'
local setups = require 'setups'

local fieldLoopTest = {
	run_name = 'FieldLoopTest',
	output_directory = 'data/field_loop_test',
	final_time = 2.0,
	checkpoint_interval = 1.0,
	vtk_output_interval = 0.01,
	cfl_parameter = 0.3,
	initial_data = setups['mhd_field_loop'],
	grid_geometry = 'cartesian',
	resolution = {48, 48, 1},
	domain_lower = {-0.5, -0.5, -0.5},
	domain_upper = { 0.5,  0.5,  0.5},
	conservation_law = {'newtonian_mhd'},
	riemann_solver = 'hlle',
	flux_scheme = {'method_of_lines_plm', plm_theta=1.5},
	runge_kutta_order = 2,
	boundary_condition = 'periodic',
}

mara.run (fieldLoopTest)
