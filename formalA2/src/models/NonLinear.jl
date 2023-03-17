
function NLtwoVarSys1(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        mulsub(0.01, q[2], mulT(q[1], q[1], cache[2]), cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(negateT(q[1], cache[2]), q[2], 2.0, cache[1])
        #= none:1 =#
        return nothing
    end
end
function lorenz(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        mulT(10.0, subT(q[2], q[1], cache[2]), cache[1])
        #= none:1 =#
        return nothing
    elseif #= none:1 =# j == 2
        #= none:1 =#
        mulsub(q[1], subT(28.0, q[3], cache[2]), q[2], cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        mulsub(q[1], q[2], mulT(2.66667, q[3], cache[2]), cache[1])
        #= none:1 =#
        return nothing
    end
end

function loktavoltera(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        mulsub(2.0, q[1], mulTT(1.2, q[1], q[2], cache[2], cache[3]), cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        muladdT(-3.0, q[2], mulT(q[1], q[2], cache[2]), cache[1])
        #= none:1 =#
        return nothing
    end
end
function vanderpol(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        createT(q[2], cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        mulsub(subT(1.0, mulT(q[1], q[1], cache[3]), cache[2]), q[2], q[1], cache[1])
        #= none:1 =#
        return nothing
    end
end
function oregonator(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        mulT(100.8, subadd(mulT(9.523809523809524e-5, q[2], cache[3]), mulT(q[1], q[2], cache[4]), mulT(q[1], subT(1.0, q[1], cache[6]), cache[5]), cache[2]), cache[1])
        #= none:1 =#
        return nothing
    elseif #= none:1 =# j == 2
        #= none:1 =#
        mulT(40320.0, subadd(mulT(-9.523809523809524e-5, q[2], cache[3]), mulT(q[1], q[2], cache[4]), q[3], cache[2]), cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subT(q[1], q[3], cache[1])
        #= none:1 =#
        return nothing
    end
end
