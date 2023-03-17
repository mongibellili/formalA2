

function foo()
    du[1] = u[2] - u[1]
    du[2] = u[3] - u[2]
    du[3] = u[4] - u[3]
    du[4] = u[5] - u[4]
    du[5] = u[6] - u[5]
end
function foo()
    du[1] = 3 * 1 * u[2] - 2 * u[1]
    du[2] = 3 * 2 * u[3] - 2 * u[2]
    du[3] = 3 * 3 * u[4] - 2 * u[3]
    du[4] = 3 * 4 * u[5] - 2 * u[4]
    du[5] = 3 * 5 * u[6] - 2 * u[5]
end
function foo()
    du[1] = ((U / Rs[1] - iL[1]) * ((Rs[1] * Rd[1]) / (Rs[1] + Rd[1])) - uC) / L
    du[2] = ((U / Rs[2] - iL[2]) * ((Rs[2] * Rd[2]) / (Rs[2] + Rd[2])) - uC) / L
    du[3] = ((U / Rs[3] - iL[3]) * ((Rs[3] * Rd[3]) / (Rs[3] + Rd[3])) - uC) / L
    du[4] = ((U / Rs[4] - iL[4]) * ((Rs[4] * Rd[4]) / (Rs[4] + Rd[4])) - uC) / L
end
relaxedqssA.NLODEProblem{2, 1, 0, 0}(3, [1.0, 1.0], [0.0], StaticArrays.SVector{2, SymEngine.Basic}[[2.0 - 1.2*q2, -1.2*q1], [q2, -3.0 + q1]], :(function f(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
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
  end), Expr[], Expr[], StaticArrays.SVector{1, SymEngine.Basic}[[0], [0]], StaticArrays.SVector{2, SymEngine.Basic}[], StaticArrays.SVector{1, SymEngine.Basic}[], relaxedqssA.EventDependencyStruct[])
