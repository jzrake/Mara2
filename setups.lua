local setups = {}


-- One profile for each component of the scalar advection law.
function setups.scalar_pulse(x, y, z)
	return {math.exp(-(x^2 + y^2) / 0.025), 2 - x*x}
end


-- Uniform hydrodynamic state with no fluid motion.
function setups.hydro_uniform(x, y, z)
	return {1.0, 0.0, 0.0, 0.0, 1.0}
end


-- A gaussian density enhancement that advects to the right with speed 1.
function setups.hydro_pulse(x, y, z)
	local d = 1.0 + math.exp(-(x^2 + y^2 + z^2) / 0.025)
	local u = 1.0
	local v = 0.0
	local w = 0.0
	local p = 1.0
	return {d, u, v, w, p}
end


-- Classic hydro shocktube setup.
function setups.shocktube1(x, y, z)
	local d = x < 0.0 and 1.0 or 0.1
	local u = 0.0
	local v = 0.0
	local w = 0.0
	local p = x < 0.0 and 1.0 or 0.125
	return {d, u, v, w, p}
end


-- Classic MHD shocktube setup.
function setups.shocktube1_mhd(x, y, z)
	local d = x < 0.0 and 1.0 or 0.1
	local u = 0.0
	local v = 0.0
	local w = 0.0
	local p = x < 0.0 and 1.0 or 0.125
	local b1 = 1.0
	local b2 = 0.0
	local b3 = 0.0
	return {d, u, v, w, p, b1, b2, b3}
end


-- Left and right states have the same total pressure, so the solution should
-- remain near the initial value.
function setups.mhd_pressure_equilibrium(x, y, z)
	local d = 1.0
	local u = 0.0
	local v = 0.0
	local w = 0.0
	local p = x < 0.0 and 1.0 or 0.5
	local b1 = 0.0
	local b2 = x < 0.0 and 0.0 or 1.0
	local b3 = 0.0
	return {d, u, v, w, p, b1, b2, b3}
end


-- Unform vertical magnetic field.
function setups.mhd_uniform_vertical(x, y, z)
	local vx = 0.0 * z
	local bz = 0.5
	return {1.0, vx, 0.0, 0.0, 1.0, 0.0, 0.0, bz}
end


-- Force-free equilibrium state.
function setups.mhd_abc_equilibrium(x, y, z)
	local A = 1.0
	local B = 1.0
	local C = 1.0
	local k = 4.0 * math.pi
	local d = 1.0
	local u = 0.0
	local v = 0.0
	local w = 0.0
	local p = 1.0
	local b1 = C * math.cos(k * z) - B * math.sin(k * y)
	local b2 = A * math.cos(k * x) - C * math.sin(k * z)
	local b3 = B * math.cos(k * y) - A * math.sin(k * x)
	return {d, u, v, w, p, b1, b2, b3}
end


return setups
