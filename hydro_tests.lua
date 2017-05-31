local mara = require 'mara'

local densityWaveNumber = 5.0
local densityAmplitude = 5e-1

local function density_loop(x, y, z)
	local R = (x^2 + y^2)^0.5
	local d = 1.0 + densityAmplitude * math.exp (-(densityWaveNumber * R)^3)
	local p = 1.0
	return {d, 1.0, 0.0, 0.0, p}
end

local function density_wave(x, y, z)
	local d = 1.0 + 0.5 * math.sin(densityWaveNumber * math.pi * x)
	local p = 1.0
	return {d, 1.0, 0.0, 0.0, p}
end

local setup = {
	run_name = 'DensityWaveTest',
	output_directory = 'data/density_wave',
	initial_data = density_wave,
	final_time = 0.0,
	checkpoint_interval = 1.0,
	vtk_output_interval = 0,--0.025,
	vtk_use_binary = true,
	cfl_parameter = 0.3,
	grid_geometry = 'cartesian',
	resolution = {99, 1, 1},
	domain_lower = {-0.5, -0.5, -0.5},
	domain_upper = { 0.5,  0.5,  0.5},
	conservation_law = {'newtonian_hydro'},
	riemann_solver = 'hlle',
	flux_scheme = {'method_of_lines_plm', plm_theta=1.5},
	runge_kutta_order = 2,
	boundary_condition = 'periodic',
}

mara.run(setup)
