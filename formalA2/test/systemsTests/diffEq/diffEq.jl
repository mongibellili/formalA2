




using OrdinaryDiffEq
#using qssv01

using Plots
function odeDiffEquPackage()
    function f(du,u,p,t)
       
        
        du[1] = -0.75*u[1]+0.2*u[2]+1.0
                du[2] =-12.65*u[1]+2.0*u[2]+1.0
       
    end
   
    u0  = [-1.0, -2.0]
    tspan = (0.0,3.0)
    #p = -9.8
    prob = ODEProblem(f,u0,tspan)
    #sol = solve(prob,Tsit5(),callback=cb)
    sol = solve(prob,BS3())
    p1=plot!(sol,marker=(:circle),markersize=2)
   # p1=plot!(sol,marker=(:circle),markersize=2,xlims=(0.0,30.0) ,ylims=(-2.04e-1,2.0))
  # p1=plot!(sol,marker=(:circle),markersize=2,xlims=(0.0,30.0) ,ylims=(-2.04e-1,2.06e-1))
   savefig(p1, "bs3_K1_ft3")
 end
#@btime
 odeDiffEquPackage() 