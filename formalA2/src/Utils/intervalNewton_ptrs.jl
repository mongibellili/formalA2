#using BenchmarkTools

@inline function smallest(args::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds m = args[1]
	for i in 2:N
		@inbounds arg = args[i]
		m = arg < m ? arg : m
	end
	m
end

@inline function horner(x::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds v = coeffs[1]
	for i in 2:N
		@inbounds coeff = coeffs[i]
		v = muladd(x, v, coeff)
	end
	v
end

@inline function horner2(x::F, y::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds v = coeffs[1]
	w = v
	for i in 2:N
		@inbounds coeff = coeffs[i]
		v = muladd(x, v, coeff)
		w = muladd(y, w, coeff)
	end
	v, w
end

@inline function hornerd(x::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds v = coeffs[1]
	d = zero(F)
	for i in 2:N
		@inbounds coeff = coeffs[i]
		d = muladd(x, d, v)
		v = muladd(x, v, coeff)
	end
	v, d
end

@inline function intervalhorner(low::F, high::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds colow = coeffs[1]
	cohigh = colow
	if low < zero(F)
		for i in 2:N
			@inbounds coeff = coeffs[i]
			colow, cohigh = if colow > zero(F)
				muladd(low, cohigh, coeff), muladd(high, colow, coeff)
			elseif cohigh < zero(F)
				muladd(high, cohigh, coeff), muladd(low, colow, coeff)
			else
				muladd(low, cohigh, coeff), muladd(low, colow, coeff)
			end
		end
	else
		for i in 2:N
			@inbounds coeff = coeffs[i]
			colow, cohigh = if colow > zero(F)
				muladd(low, colow, coeff), muladd(high, cohigh, coeff)
			elseif cohigh < zero(F)
				muladd(high, colow, coeff), muladd(low, cohigh, coeff)
			else
				muladd(high, colow, coeff), muladd(high, cohigh, coeff)
			end
		end	
	end
    colow, cohigh
end

@inline function posintervalhorner(low::F, high::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds colow = coeffs[1]
	cohigh = colow
	for i in 2:N
		@inbounds coeff = coeffs[i]
		colow, cohigh = if colow > zero(F)
			muladd(low, colow, coeff), muladd(high, cohigh, coeff)
		elseif cohigh < zero(F)
			muladd(high, colow, coeff), muladd(low, cohigh, coeff)
		else
			muladd(high, colow, coeff), muladd(high, cohigh, coeff)
		end
	end
	colow, cohigh
end

function smallestpositiverootintervalnewtonregulafalsi(coeffs::NTuple{N,F}, doms::Ptr{NTuple{2,F}}) where {N, F <: AbstractFloat}
	if N == 1
		return typemax(F)
	elseif N == 2
		@inbounds ret = -coeffs[2] / coeffs[1]
		if ret < zero(F)
			return typemax(F)
		else
			return ret
		end
	end
    @inbounds _coeff1 = inv(coeffs[1])
    poly = ntuple(N) do i
        @inbounds _coeff1 * coeffs[i]
    end
    poly′ = ntuple(N - 1) do i
        @inbounds (N-i) * poly[i]
    end
	MM = smallest(poly)
	if MM > zero(F) return typemax(F) end
	domlow, domhigh = zero(F), one(F) - MM
	index = 0
    while true
		#@show domlow domhigh
		mid = 0.5(domlow + domhigh)
		comid = horner(mid, poly)
		codom′low, codom′high = posintervalhorner(domlow, domhigh, poly′)
		#@show comid codom′low codom′high
		if codom′low < zero(F) < codom′high
			leftlow, lefthigh, rightlow, righthigh = if comid < zero(F)
				domlow, mid - comid / codom′low, mid - comid / codom′high, domhigh
			else
				domlow, mid - comid / codom′high, mid - comid / codom′low, domhigh
			end
			#@show leftlow lefthigh rightlow righthigh
			if leftlow < lefthigh
				if rightlow < righthigh
					index += 1
					unsafe_store!(doms, (rightlow, righthigh), index)
				end
				domlow, domhigh = leftlow, lefthigh
				continue
			elseif rightlow < righthigh
				domlow, domhigh = rightlow, righthigh
				continue
			end
		else
			codomlow, codomhigh = horner2(domlow, domhigh, poly)
			#@show domlow domhigh codomlow codomhigh
			if codomlow * codomhigh < zero(F)
				while true
					comid, comid′ = hornerd(mid, poly)
					delta = comid / comid′
					newmid = mid - delta
					if abs(delta) < 1.0e-6mid 
						return newmid
					elseif domlow < newmid < domhigh
						mid = newmid
					else
						if comid * codomlow < zero(F)
							domhigh, codomhigh = mid, comid
						else
							domlow, codomlow = mid, comid
						end
						mid = (domlow*codomhigh - domhigh*codomlow) / (codomhigh - codomlow) # regula falsi
					end
				end
			end
		end
		if index == 0 break end
		domlow, domhigh = unsafe_load(doms, index)
		index -= 1
	end
	return typemax(F)
end
#coeffs=NTuple{4,Float64}((-1.33,2.0,1.0,6.0))
#= coeffs=NTuple{4,Float64}((1.33,12.65,1.0,-6.0))
[0.00E+00	-0.2906652780025	-8938.3815305787700	-51442410.0831105000000] =#
#= pp=pointer(Vector{NTuple{2,Float64}}(undef, 3))
#pp=pointer(Vector{NTuple{2,Float64}}(undef, 1))#signal (11): Segmentation fault
#@show pp
#@btime gensmallestpositiverootintervalnewton(1.0, -2.0, -1.0)#64.912 ns (1 allocation: 32 bytes)
#@show gensmallestpositiverootintervalnewton(1.0, -2.0, -1.0)#2.414213562506035#@show quadraticsmallestpositiveroot(1.0, -2.0, -1.0)#1.649 ns (0 allocations: 0 bytes)
#@show smallestpositiverootintervalnewtonregulafalsi((0.00E+00	,-0.2906652780025	,-8938.3815305787700	,-51442410.0831105000000),pp)
d=1e-6
#c= 0.4289402307079971-0.46926407070799725
c=-0.040323840000000166
#b=7597.313022054603-11511.429846502817
b=-3914.116824448214
a=1.021243315770821e8
@show smallestpositiverootintervalnewtonregulafalsi((a, b, c, 1e-6),pp) =#
#@show smallestpositiverootintervalnewtonregulafalsi((1e-6, c, b, a),pp)
#@show smallestpositiverootintervalnewtonregulafalsi((0.0, -0.2906652780025, -8938.3815305787700, -5.144241008311048e7),pp)
#@show smallestpositiverootintervalnewtonregulafalsi((-1.33,2.0,1.0,6.0),pp)
#@btime smallestpositiverootintervalnewtonregulafalsi((-1.33,2.0,1.0,6.0),pp)# 93.687 ns (0 allocations: 0 bytes)
#@btime smallestpositiverootintervalnewtonregulafalsi((-1.33,2.0,1.0,6.0),pp)#126.486 ns (2 allocations: 64 bytes)....2.5155562382541037...cubic:18ns
#@btime gensmallestpositiverootintervalnewton(1.33,2.0,-2.63,6.0)#32.639 ns (0 allocations: 0 bytes) ....Inf
#@show gensmallestpositiverootintervalnewton(2.0,-3.0,-0.5,-5.0)
#@btime gensmallestpositiverootintervalnewton(-1.33,2.0,1.0,6.0,2.2)#171.461 ns (2 allocations: 64 bytes)
#@btime gensmallestpositiverootintervalnewton(-1.33,2.0,1.0,6.0,2.2,3.14)#220.877 ns (2 allocations: 64 bytes)
