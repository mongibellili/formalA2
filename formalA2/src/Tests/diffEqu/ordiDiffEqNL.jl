#using OrdinaryDiffEq
using DifferentialEquations
#using qssv01
using BenchmarkTools
using Plots
function odeDiffEquPackage()
 
    function funcName(du,u,p,t)# api requires four args

    
###############lorenz
        #= du[1] = 10.0*(u[2]-u[1])
        du[2] = u[1]*(28.0-u[3]) - u[2]
        du[3] = u[1]*u[2] - (2.66667)*u[3] =#
##################lokta voltera
       #=  du[1] = 2.0 * u[1] - 1.2 * u[1]*u[2]
      du[2] = -3.0 * u[2] + u[1]*u[2] =#
##############Van der pol
      #= du[1] = u[2]
      du[2] = (1-u[1]*u[1]) * u[2] - u[1] =#
###################oregonator
#= eta1_default = 0.00992063492063492
eta2_default = 2.48015873015873e-5
q_default = 9.523809523809524e-5 =#




du[1] = 100.8*(   9.523809523809524e-5 * u[2] - u[1] * u[2] + u[1] * ( 1.0 - u[1] ) ) 
du[2] = 40320.0*( - 9.523809523809524e-5 * u[2] - u[1] * u[2] +  u[3] ) 
du[3] = u[1] - u[3]

#= du[1] =    0.0096*u[2] - 100.8*u[1]*u[2] + 100.8*u[1] -100.8*u[1]* u[1] 
du[2] =  3.84*u[2] - 40320.0*u[1]*u[2] + 40320.0*u[3] 
du[3] = u[1] - u[3] =#



     #=  du[1] = 100.0*(u[1]+u[2]-u[1]*u[2]-1e-6*u[1]*u[1])
      du[2] = 0.01*(-u[2]-u[1]*u[2]+1.1*u[3]) 
      du[3] = 3.835*(u[1]-u[3]) =#

    end
    tspan = (0.0,1.0)
    u0 = [1.0,0.0,0.0]
    prob = ODEProblem(funcName,u0,tspan)
    sol = solve(prob,BS3(),saveat=0.05,abstol = 1e-6, reltol = 1e-3)
    #sol = solve(prob,Rodas5P(),saveat=1.0,abstol = 1e-6, reltol = 1e-3)
    #sol = solve(prob,Rosenbrock23(),abstol = 1e-6, reltol = 1e-3)
  
   
   
   


   # p1=plot(sol,idxs = (0, 1)#= ,marker=(:circle),ylims=(0,1) =#)
    display(plot(sol, idxs = (0, 2),legend = false))
    println("press enter to exit")
  readline()
   # p1=plot!(sol[2],marker=(:circle),ylims=(0,1))
   # display(plot!(sol,xlims=(100,160),ylims=(-0.000002,0.000002)))
  #  display(plot!(sol,xlims=(10.0,10.0005),ylims=(-0.49971,-0.49960)))
  #savefig(p1, "plot_sysOREGONATOR_bs3_1e-3_ft20.png")
    
end
#@time 
odeDiffEquPackage()  

