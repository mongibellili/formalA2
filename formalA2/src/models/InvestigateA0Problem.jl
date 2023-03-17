function A011(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        mulT(2.0, q[2], cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        mulsub(-0.1, q[1], mulT(10.0, q[2], cache[2]), cache[1])
        #= none:1 =#
        return nothing
    end
end
function A012(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        mulsub(2.0, q[2], 2.0, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(mulT(-0.1, q[1], cache[2]), mulT(10.0, q[2], cache[3]), 10.0, cache[1])
        #= none:1 =#
        return nothing
    end
end
function A013(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        mulsub(2.0, q[2], 40.0, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(mulT(-0.1, q[1], cache[2]), mulT(10.0, q[2], cache[3]), 200.0, cache[1])
        #= none:1 =#
        return nothing
    end
end
function A014(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        mulsub(2.0, q[2], 0.02, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(mulT(-0.1, q[1], cache[2]), mulT(10.0, q[2], cache[3]), 0.1, cache[1])
        #= none:1 =#
        return nothing
    end
end
function A015(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        mulsub(2.0, q[2], 2.0, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(mulT(-0.1, q[1], cache[2]), mulT(10.0, q[2], cache[3]), 9.8, cache[1])
        #= none:1 =#
        return nothing
    end
end
function A021(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        mulT(0.25, q[2], cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        mulsub(-8.0, q[1], mulT(4.0, q[2], cache[2]), cache[1])
        #= none:1 =#
        return nothing
    end
end
function A022(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        mulsub(0.25, q[2], 0.25, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(mulT(-8.0, q[1], cache[2]), mulT(4.0, q[2], cache[3]), 4.0, cache[1])
        #= none:1 =#
        return nothing
    end
end
function A023(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        mulsub(0.25, q[2], 5.0, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(mulT(-8.0, q[1], cache[2]), mulT(4.0, q[2], cache[3]), 80.0, cache[1])
        #= none:1 =#
        return nothing
    end
end
function A024(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        mulsub(0.25, q[2], 0.0025, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(mulT(-8.0, q[1], cache[2]), mulT(4.0, q[2], cache[3]), 0.04, cache[1])
        #= none:1 =#
        return nothing
    end
end
function A025(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        mulsub(0.25, q[2], 0.25, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subsub(mulT(-8.0, q[1], cache[2]), mulT(4.0, q[2], cache[3]), 12.0, cache[1])
        #= none:1 =#
        return nothing
    end
end
