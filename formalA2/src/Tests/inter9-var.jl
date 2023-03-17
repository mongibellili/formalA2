using BenchmarkTools
using StaticArrays

@inline function smallest(args::NTuple{N,Float64}) where { N }
	@inbounds m = args[1]
	for i in 2:N
		@inbounds arg = args[i]
		m = arg < m ? arg : m
	end
	m
end

@inline function horner(x::Float64, coeffs::NTuple{N,Float64}) where { N }
	@inbounds v = coeffs[1]
	for i in 2:N
		@inbounds coeff = coeffs[i]
		v = muladd(x, v, coeff)
	end
	v
end

@inline function horner2(x::Float64, y::Float64, coeffs::NTuple{N,Float64}) where { N }
	@inbounds v = coeffs[1]
	w = v
	for i in 2:N
		@inbounds coeff = coeffs[i]
		v = muladd(x, v, coeff)
		w = muladd(y, w, coeff)
	end
	v, w
end

@inline function hornerd(x::Float64, coeffs::NTuple{N,Float64}) where { N }
	@inbounds v = coeffs[1]
	d =  0.0
	for i in 2:N
		@inbounds coeff = coeffs[i]
		d = muladd(x, d, v)
		v = muladd(x, v, coeff)
	end
	v, d
end

@inline function intervalhorner(low::Float64, high::Float64, coeffs::NTuple{N,Float64}) where { N }
	@inbounds colow = coeffs[1]
	cohigh = colow
	if low <  0.0
		for i in 2:N
			@inbounds coeff = coeffs[i]
			colow, cohigh = if colow >  0.0
				muladd(low, cohigh, coeff), muladd(high, colow, coeff)
			elseif cohigh <  0.0
				muladd(high, cohigh, coeff), muladd(low, colow, coeff)
			else
				muladd(low, cohigh, coeff), muladd(low, colow, coeff)
			end
		end
	else
		for i in 2:N
			@inbounds coeff = coeffs[i]
			colow, cohigh = if colow >  0.0
				muladd(low, colow, coeff), muladd(high, cohigh, coeff)
			elseif cohigh <  0.0
				muladd(high, colow, coeff), muladd(low, cohigh, coeff)
			else
				muladd(high, colow, coeff), muladd(high, cohigh, coeff)
			end
		end	
	end
    colow, cohigh
end

@inline function posintervalhorner(low::Float64, high::Float64, coeffs::NTuple{N,Float64}) where { N }
	@inbounds colow = coeffs[1]
	cohigh = colow
	for i in 2:N
		@inbounds coeff = coeffs[i]
		colow, cohigh = if colow >  0.0
			muladd(low, colow, coeff), muladd(high, cohigh, coeff)
		elseif cohigh <  0.0
			muladd(high, colow, coeff), muladd(low, cohigh, coeff)
		else
			muladd(high, colow, coeff), muladd(high, cohigh, coeff)
		end
	end
	colow, cohigh
end

function smallestpositiverootintervalnewtonregulafalsi(coeffs::NTuple{N,Float64}, doms::Ptr{NTuple{2,Float64}}) where { N }
	if N == 1
		return Inf
	elseif N == 2
		@inbounds ret = -coeffs[2] / coeffs[1]
		if ret <  0.0
			return Inf
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
	if MM >  0.0 return Inf end
	domlow, domhigh =  0.0, one(Float64) - MM
	index = 0
    while true
		#@show domlow domhigh
		mid = 0.5(domlow + domhigh)
		comid = horner(mid, poly)
		codom′low, codom′high = posintervalhorner(domlow, domhigh, poly′)
		#@show comid codom′low codom′high
		if codom′low <  0.0 < codom′high
			leftlow, lefthigh, rightlow, righthigh = if comid <  0.0
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
			if codomlow * codomhigh <  0.0
				while true
					comid, comid′ = hornerd(mid, poly)
					delta = comid / comid′
					newmid = mid - delta
					if abs(delta) < 1.0e-2mid 
						return newmid
					elseif domlow < newmid < domhigh
						mid = newmid
					else
						if comid * codomlow <  0.0
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
	return Inf
end



function minPosRoot(coeff::SVector{2,Float64}, ::Val{1}) # coming from val(1) means coef has x and derx only...size 2
    mpr=-1
        if coeff[2] == 0 
            mpr = Inf
        else 
            mpr = -coeff[1] / coeff[2];
        end
        if mpr < 0
            mpr = Inf
        end
       # println("mpr inside minPosRoot in utils= ",mpr)
    return mpr
end
function minPosRoot(coeff::SVector{3,Float64}, ::Val{2}) # credit goes to github.com/CIFASIS/qss-solver
    mpr=-1 #coef1=c, coef2=b, coef3=a
    if coeff[3] == 0 #|| (10000 * abs(coeff[3])) < abs(coeff[2])
        if coeff[2] == 0
          mpr = Inf
        else 
          mpr = -coeff[1] / coeff[2]
        end
        if mpr < 0
          mpr = Inf
        end
    else  #double disc;
        disc = coeff[2] * coeff[2] - 4 * coeff[3] * coeff[1]#b^2-4ac
        if disc < 0 # no real roots
          mpr = Inf
        else  #double sd, r1;
          sd = sqrt(disc);
          r1 = (-coeff[2] + sd) / (2 * coeff[3]);
          if r1 > 0 
            mpr = r1;
          else 
            mpr = Inf;
          end
          r1 = (-coeff[2] - sd) / (2 * coeff[3]);
          if ((r1 > 0) && (r1 < mpr)) 
            mpr = r1;
          end
        end       
    end
    return mpr
end


@inline function minPosRoot(coeffs::SVector{4,Float64}, ::Val{3})#where F <: AbstractFloat
  if coeffs[4] == 0.0
    coeffs2=@SVector[coeffs[1],coeffs[2],coeffs[3]]
    return minPosRoot(coeffs2, Val(2))
  end
  _a = 1.0 / coeffs[4]
  b, c, d = coeffs[3] * _a, coeffs[2] * _a, coeffs[1] * _a
  m = b < c ? b : c
  m = d < m ? d : m
  m > eps(Float64) && return Inf#typemax(Float64) # Cauchy bound
  _3 = 1.0 / 3
  _9 = 1.0 / 9
  SQ3 = sqrt(3.0)
  xₙ = -b * _3
  b²_9 = b * b * _9
  yₙ = muladd(muladd(-2, b²_9, c), xₙ, d)   #eq to 2R
  δ² = muladd(-_3, c, b²_9)                  #eq to Q
  h² = 4δ² * δ² * δ²
  Δ = muladd(yₙ, yₙ, -h²)
  if Δ > 4eps(Float64) # one real root and two complex roots
  p = yₙ < 0 ? cbrt(0.5 * (-yₙ + √Δ)) : cbrt(0.5 * (-yₙ - √Δ))
  q = δ² / p
  z = xₙ + p + q
  z > -eps(Float64) ? z : Inf#typemax(Float64)
  elseif Δ < -4eps(Float64) # three real roots
  θ = abs(yₙ) < eps(Float64) ? 0.5π * _3 : atan(√abs(Δ) / abs(yₙ)) * _3 # acos(-yₙ / √h²)
  δ = yₙ < 0 ? √abs(δ²) : -√abs(δ²)
  z₁ = 2δ * cos(θ)
  z₂ = muladd(-0.5, z₁, xₙ)
  z₃ = SQ3 * δ * sin(θ)
  x₁ = xₙ + z₁
  x₂ = z₂ + z₃
  x₃ = z₂ - z₃
  x = x₁ > -eps(Float64) ? x₁ : Inf# typemax(F)
  x = x₂ > -eps(Float64) && x₂ < x ? x₂ : x
  x₃ > -eps(Float64) && x₃ < x ? x₃ : x
  else # double or triple real roots
  δ = cbrt(0.5yₙ)
  x₁ = xₙ + δ
  x₂ = xₙ - 2δ
  x = x₁ > -eps(Float64) ? x₁ : Inf#typemax(F)
  x₂ > -eps(Float64) && x₂ < x ? x₂ : x
  end
end

pp=pointer(Vector{NTuple{2,Float64}}(undef, 1))
#-----------------quadratic-----------------------------
coeffsQuad=NTuple{3,Float64}((12.65,1.0,-6.0))
 #=@show smallestpositiverootintervalnewtonregulafalsi(coeffsQuad,pp)
@btime smallestpositiverootintervalnewtonregulafalsi(coeffsQuad,pp) =#
#=@show  smallestpositiverootintervalnewtonregulafalsi((12.65,1.0,-6.0),pp)
@time  smallestpositiverootintervalnewtonregulafalsi((12.65,1.0,-6.0),pp)
@timev smallestpositiverootintervalnewtonregulafalsi((12.65,1.0,-6.0),pp)
@btime smallestpositiverootintervalnewtonregulafalsi((12.65,1.0,-6.0),pp)
 coeffs=@SVector[-6.0,1.0,12.65]
@show minPosRoot(coeffs,Val(2))
@time minPosRoot(coeffs,Val(2))
@timev minPosRoot(coeffs,Val(2))
@btime minPosRoot(coeffs,Val(2)) =#
#-----------------cubic-----------------------------
#= coeffsCubic=NTuple{4,Float64}((1.33,12.65,1.0,-6.0))
 @show smallestpositiverootintervalnewtonregulafalsi(coeffs,pp)
@btime smallestpositiverootintervalnewtonregulafalsi(coeffs,pp) =#
@show smallestpositiverootintervalnewtonregulafalsi((1.33,12.65,1.0,-6.0),pp)
 @time smallestpositiverootintervalnewtonregulafalsi((1.33,12.65,1.0,-6.0),pp)
@timev smallestpositiverootintervalnewtonregulafalsi((1.33,12.65,1.0,-6.0),pp)
@btime smallestpositiverootintervalnewtonregulafalsi((1.33,12.65,1.0,-6.0),pp)
coeffs=@SVector[-6.0,1.0,12.65,1.33]
@show minPosRoot(coeffs,Val(3))
@time minPosRoot(coeffs,Val(3)) 
@timev minPosRoot(coeffs,Val(3)) 
@btime minPosRoot(coeffs,Val(3)) 