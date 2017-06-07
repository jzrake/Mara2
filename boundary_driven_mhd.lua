local vAlfven = 1.0
local vDrive  = 0.2


math.randomseed(2)      -- The RNG may be used to generate noise
local noise_level = 0.0 -- May be overwritten in apply_model

local function background(x, y, z)
    local d = 1.0 + noise_level * (math.random() - 0.5)
    return {d, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, vAlfven}
end

local function abc_driving_floor_ceil(x, y, z)
    local sgnz = z > 0 and 1 or -1
    local vx = vDrive * math.sin(4 * math.pi * y) * sgnz
    local vy = vDrive * math.cos(4 * math.pi * x) * sgnz
    return {vx, vy}
end

local function abc_driving_floor_only(x, y, z)
    if (z < 0) then
        local vx = vDrive * math.sin(4 * math.pi * y)
        local vy = vDrive * math.cos(4 * math.pi * x)
        return {vx, vy}
    else
        return {0, 0}
    end
end

local function single_vortex_pair(x, y, z)
    local sgnz = z > 0 and 1 or -1
    local R = (x^2 + y^2)^0.5
    local k = 10
    local v = vDrive / k
    local vf = v * 3 * k * (k * R)^2 * math.exp (-(k * R)^3) -- v-phi
    local vx = vf * (-y / R) * sgnz
    local vy = vf * ( x / R) * sgnz
    return {vx, vy}
end

local function apply_model(id, res)

    local models = {
        ['bdrv-noceilz2'] =
        {
            bvf         = abc_driving_floor_only,
            aspect      = 2
        },
        ['bdrv-noceilz4'] =
        {
            bvf         = abc_driving_floor_only,
            aspect      = 4
        },
        ['bdrv-noceilz2-noise2'] =
        {
            bvf         = abc_driving_floor_only,
            aspect      = 2,
            noise_level = 1e-2
        },
    }

    local M = models[id]
    noise_level       = M.noise_level -- local variable
    run_name          = id .. '-' .. tostring(res)
    output_directory  = 'data/' .. run_name
    boundary_value    = M.bvf
    resolution        = { res, res, res * M.aspect }
    domain_lower      = {-0.5, -0.5, -0.5 * M.aspect }
    domain_upper      = { 0.5,  0.5,  0.5 * M.aspect }
end

final_time            = 72
initial_data          = background
checkpoint_interval   = 0.50
vtk_output_interval   = 0.00
time_series_interval  = 0.01
cfl_parameter         = 0.33
mesh_geometry         = 'cartesian'
conservation_law      = 'newtonian_mhd'
riemann_solver        = 'hlle'
flux_scheme           = 'method_of_lines_plm'
plm_theta             = 1.75
runge_kutta_order     = 2
boundary_condition    = 'driven_mhd'
disable_ct            = false
constrained_transport = 'uniform_cartesian'
pressure_floor        = 1e-2

apply_model('bdrv-noceilz2-noise2', 32)
