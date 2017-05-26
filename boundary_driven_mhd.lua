local mara = require 'mara'

local fieldMagnitude = 1.0

local function background(x, y, z)
	return {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, fieldMagnitude}
end

local setup = {
	run_name = 'BoundaryDrivenMHD',
	output_directory = 'data/bdrv-mhd-16',
	final_time = 64.,
	checkpoint_interval = 0.0,
	vtk_output_interval = 0.1,
	cfl_parameter = 0.3,
	initial_data = background,
	grid_geometry = 'cartesian',
	resolution = {16, 16, 16},
	domain_lower = {-0.5, -0.5, -0.5},
	domain_upper = { 0.5,  0.5,  0.5},
	conservation_law = {'newtonian_mhd'},
	riemann_solver = 'hlle',
	flux_scheme = {'method_of_lines_plm', plm_theta=1.5},
	runge_kutta_order = 3,
	boundary_condition = 'driven_mhd',
	disable_ct = false,
	pressure_floor = 1e-2
}

mara.run(setup)
