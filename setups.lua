local setups = {}

function setups.scalar_pulse(x, y, z)
	return {math.exp(-(x^2 + y^2) / 0.025), 2 - x*x}
end

function setups.euler_pulse(x, y, z)
	local d = 1.0 + math.exp(-(x^2 + y^2 + z^2) / 0.025)
	local u = 0.0
	local v = 0.0
	local w = 0.0
	local p = 1.0
	return {d, u, v, w, p}
end

function setups.shocktube1(x, y, z)
	local d = x^2 + y^2 + z^2 < 0.1 and 1.0 or 0.1
	local u = 0.0
	local v = 0.0
	local w = 0.0
	local p = x^2 + y^2 + z^2 < 0.1 and 1.0 or 0.125
	return {d, u, v, w, p}
end

return setups
