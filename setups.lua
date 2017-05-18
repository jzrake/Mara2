local setups = {}

function setups.scalar_pulse(x, y, z)
	return {math.exp(-(x^2 + y^2) / 0.025), 2 - x*x}
end

function setups.hydro_pulse(x, y, z)
	local d = 1.0 + math.exp(-(x^2 + y^2 + z^2) / 0.025)
	local u = 1.0
	local v = 0.0
	local w = 0.0
	local p = 1.0
	return {d, u, v, w, p}
end

function setups.shocktube1(x, y, z)
	local d = x < 0.0 and 1.0 or 0.1
	local u = 0.0
	local v = 0.0
	local w = 0.0
	local p = x < 0.0 and 1.0 or 0.125
	return {d, u, v, w, p}
end

function setups.shocktube1_mhd(x, y, z)
	local d = x < 0.0 and 1.0 or 0.1
	local u = 0.0
	local v = 0.0
	local w = 0.0
	local p = x < 0.0 and 1.0 or 0.125
	local b1 = 0.0
	local b2 = 0.0
	local b3 = 0.0
	return {d, u, v, w, p, b1, b2, b3}
end


return setups
