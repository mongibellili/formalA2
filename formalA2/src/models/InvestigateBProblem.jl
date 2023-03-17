function B11(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        muladdT(0.15, q[1], mulT(0.25, q[2], cache[2]), cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        mulsub(-7.0, q[1], mulT(4.0, q[2], cache[2]), cache[1])
        #= none:1 =#
        return nothing
    end
end
function B12(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        addsub(mulT(0.15, q[1], cache[2]), mulT(0.25, q[2], cache[3]), 0.25, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(mulT(-7.0, q[1], cache[2]), mulT(4.0, q[2], cache[3]), 4.0, cache[1])
        #= none:1 =#
        return nothing
    end
end
function B13(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        addsub(mulT(0.15, q[1], cache[2]), mulT(0.25, q[2], cache[3]), 5.0, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(mulT(-7.0, q[1], cache[2]), mulT(4.0, q[2], cache[3]), 80.0, cache[1])
        #= none:1 =#
        return nothing
    end
end
function B14(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        addsub(mulT(0.15, q[1], cache[2]), mulT(0.25, q[2], cache[3]), 0.0025, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(mulT(-7.0, q[1], cache[2]), mulT(4.0, q[2], cache[3]), 0.04, cache[1])
        #= none:1 =#
        return nothing
    end
end
function B15(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        addT(mulT(0.15, q[1], cache[2]), mulT(0.25, q[2], cache[3]), 0.05, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subsub(mulT(-7.0, q[1], cache[2]), mulT(4.0, q[2], cache[3]), 10.0, cache[1])
        #= none:1 =#
        return nothing
    end
end
function B21(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        muladdT(0.1, q[1], mulT(0.25, q[2], cache[2]), cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        mulsub(-4.0, q[1], mulT(7.0, q[2], cache[2]), cache[1])
        #= none:1 =#
        return nothing
    end
end
function B22(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
      if j == 1
        #= none:1 =#
        addsub(mulT(0.1, q[1], cache[2]), mulT(0.25, q[2], cache[3]), 0.25, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(mulT(-4.0, q[1], cache[2]), mulT(7.0, q[2], cache[3]), 7.0, cache[1])
        #= none:1 =#
        return nothing
    end
end
function B23(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        addsub(mulT(0.1, q[1], cache[2]), mulT(0.25, q[2], cache[3]), 5.0, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(mulT(-4.0, q[1], cache[2]), mulT(7.0, q[2], cache[3]), 140.0, cache[1])
        #= none:1 =#
        return nothing
    end
end
function B24(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        addsub(mulT(0.1, q[1], cache[2]), mulT(0.25, q[2], cache[3]), 0.0025, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(mulT(-4.0, q[1], cache[2]), mulT(7.0, q[2], cache[3]), 0.07, cache[1])
        #= none:1 =#
        return nothing
    end
end
function B25(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        addsub(mulT(0.1, q[1], cache[2]), mulT(0.25, q[2], cache[3]), 0.05, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subsub(mulT(-4.0, q[1], cache[2]), mulT(7.0, q[2], cache[3]), 1.0, cache[1])
        #= none:1 =#
        return nothing
    end
end
function B31(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        muladdT(-0.25, q[1], mulT(0.15, q[2], cache[2]), cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        mulsub(-4.0, q[1], mulT(7.0, q[2], cache[2]), cache[1])
        #= none:1 =#
        return nothing
    end
end
function B32(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        addsub(mulT(-0.25, q[1], cache[2]), mulT(0.15, q[2], cache[3]), 0.15, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(mulT(-4.0, q[1], cache[2]), mulT(7.0, q[2], cache[3]), 7.0, cache[1])
        #= none:1 =#
        return nothing
    end
end
function B33(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        addsub(mulT(-0.25, q[1], cache[2]), mulT(0.15, q[2], cache[3]), 3.0, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(mulT(-4.0, q[1], cache[2]), mulT(7.0, q[2], cache[3]), 140.0, cache[1])
        #= none:1 =#
        return nothing
    end
end
function B34(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        addsub(mulT(-0.25, q[1], cache[2]), mulT(0.15, q[2], cache[3]), 0.0015, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(mulT(-4.0, q[1], cache[2]), mulT(7.0, q[2], cache[3]), 0.07, cache[1])
        #= none:1 =#
        return nothing
    end
end
function B35(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        addsub(mulT(-0.25, q[1], cache[2]), mulT(0.15, q[2], cache[3]), 0.65, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subsub(mulT(-4.0, q[1], cache[2]), mulT(7.0, q[2], cache[3]), 1.0, cache[1])
        #= none:1 =#
        return nothing
    end
end

function B41(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        muladdT(-0.01, q[1], mulT(1.24, q[2], cache[2]), cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        mulsub(-80.0, q[1], mulT(20.0, q[2], cache[2]), cache[1])
        #= none:1 =#
        return nothing
    end
end
function B42(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        addsub(mulT(-0.01, q[1], cache[2]), mulT(1.24, q[2], cache[3]), 1.24, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(mulT(-80.0, q[1], cache[2]), mulT(20.0, q[2], cache[3]), 20.0, cache[1])
        #= none:1 =#
        return nothing
    end
end
function B43(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        addsub(mulT(-0.01, q[1], cache[2]), mulT(1.24, q[2], cache[3]), 24.8, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(mulT(-80.0, q[1], cache[2]), mulT(20.0, q[2], cache[3]), 400.0, cache[1])
        #= none:1 =#
        return nothing
    end
end
function B44(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        addT(mulT(-0.01, q[1], cache[2]), mulT(1.24, q[2], cache[3]), 0.0124, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(mulT(-80.0, q[1], cache[2]), mulT(20.0, q[2], cache[3]), 0.2, cache[1])
        #= none:1 =#
        return nothing
    end
end
function B45(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        addT(mulT(-0.01, q[1], cache[2]), mulT(1.24, q[2], cache[3]), -1.22, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(mulT(-80.0, q[1], cache[2]), mulT(20.0, q[2], cache[3]), -140.0, cache[1])
        #= none:1 =#
        return nothing
    end
end
function B51(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        muladdT(-0.01, q[1], mulT(1.24, q[2], cache[2]), cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        mulsub(-80.0, q[1], mulT(20.0, q[2], cache[2]), cache[1])
        #= none:1 =#
        return nothing
    end
end
function B52(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        subadd(mulT(-20.0, q[1], cache[2]), mulT(80.0, q[2], cache[3]), 80.0, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(mulT(1.24, q[1], cache[2]), mulT(0.01, q[2], cache[3]), 0.01, cache[1])
        #= none:1 =#
        return nothing
    end
end
function B53(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        subadd(mulT(-20.0, q[1], cache[2]), mulT(80.0, q[2], cache[3]), 1600.0, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(mulT(1.24, q[1], cache[2]), mulT(0.01, q[2], cache[3]), 0.2, cache[1])
        #= none:1 =#
        return nothing
    end
end
function B54(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        subadd(mulT(-20.0, q[1], cache[2]), mulT(80.0, q[2], cache[3]), 0.8, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(mulT(1.24, q[1], cache[2]), mulT(0.01, q[2], cache[3]), 0.0001, cache[1])
        #= none:1 =#
        return nothing
    end
end
function B55(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        subadd(mulT(-20.0, q[1], cache[2]), mulT(80.0, q[2], cache[3]), 40.0, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(mulT(1.24, q[1], cache[2]), mulT(0.01, q[2], cache[3]), 2.49, cache[1])
        #= none:1 =#
        return nothing
    end
end
