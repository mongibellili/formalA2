using formalA2
#using OrdinaryDiffEq
#using BenchmarkTools
#using Plots


#include("D://QS_Solver//qss//src//models//classicProblem.jl") 
#include("/home/unknown/formalqssAB/formalqssAB/src/models/classicProblem.jl")# where you saved the model
#include("/home/unknown/formalqssAB/formalqssAB/src/models/InvestigateBProblem.jl") 
#include("D://models//NonLinear.jl")
include("/home/unknown/models/NonLinear.jl")
function test()
    odeprob = @NLodeProblem begin
      ##########lorenz equation
      #= u = [1.0, 0.0,0.0]
      discrete = [0.0]
      du[1] = 10.0*(u[2]-u[1])
      du[2] = u[1]*(28.0-u[3]) - u[2]
      du[3] = u[1]*u[2] - (2.66667)*u[3] =#
      #######the Lotka-Volterra model...
      #= u = [1.0, 1.0]
      discrete = [0.0]
      du[1] = 2.0 * u[1] - 1.2 * u[1]*u[2]
      du[2] = -3.0 * u[2] + u[1]*u[2] =#

      ###########Van der pol
      #= u = [0.0, 1.7320]
      discrete = [0.0]
      du[1] = u[2]
      du[2] = (1.0-u[1]*u[1]) * u[2] - u[1] =#
      #############oregonator
      u = [1.0, 0.0,0.0]
      discrete = [0.0]
      du[1] = 100.8*(   9.523809523809524e-5 * u[2] - u[1] * u[2] + u[1] * ( 1.0 - u[1] ) ) 
      du[2] = 40320.0*( - 9.523809523809524e-5 * u[2] - u[1] * u[2] +  u[3] ) 
      du[3] = u[1] - u[3]

  end
 # solqss2=QSS_Solve_from_model(lorenz,odeprob,10.0,mliqss3(),saveat(0.1),0.0,1e-6,1e-3)
 # solqss2=QSS_Solve_from_model(loktavoltera,odeprob,10.0,qss3(),saveat(0.1),0.0,1e-6,1e-3)
 #=  solqss2=QSS_Solve_from_model(oregonator,odeprob,10.0,qss2(),saveat(0.1),0.0,1e-6,1e-3)
  save_Sol(solqss2)
  solqss3=QSS_Solve_from_model(oregonator,odeprob,10.0,qss3(),saveat(0.1),0.0,1e-6,1e-3)
  save_Sol(solqss3) =#
  solliqss2=QSS_Solve_from_model(oregonator,odeprob,20.0,liqss2(),saveat(0.1),0.0,1e-6,1e-3)
  save_Sol(solliqss2)
  save_Sol(solliqss2,"formalold",1,2,3;xlims=(0.0,20.0),ylims=(0.0,1.005))
#=   solliqss3=QSS_Solve_from_model(oregonator,odeprob,10.0,liqss3(),saveat(0.1),0.0,1e-6,1e-3)
  save_Sol(solliqss3)
  solmliqss2=QSS_Solve_from_model(oregonator,odeprob,10.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
  save_Sol(solmliqss2)
  solmliqss3=QSS_Solve_from_model(oregonator,odeprob,10.0,mliqss3(),saveat(0.1),0.0,1e-6,1e-3)
  save_Sol(solmliqss3) =#

end

#@time 
test()



#= 
function odeDiffEquPackage()
 
    function funcName(du,u,p,t)# api requires four args

    
###############lorenz
        #= du[1] = 10.0*(u[2]-u[1])
        du[2] = u[1]*(28.0-u[3]) - u[2]
        du[3] = u[1]*u[2] - (2.66667)*u[3] =#
##################lokta voltera
        #= du[1] = 2.0 * u[1] - 1.2 * u[1]*u[2]
      du[2] = -3.0 * u[2] + u[1]*u[2] =#
##############Van der pol
      du[1] = u[2]
      du[2] = (1-u[1]*u[1]) * u[2] - u[1]

    end
    tspan = (0.0,20.0)
    u0 = [0.0, 1.7320]
    prob = ODEProblem(funcName,u0,tspan)
   # sol = solve(prob,BS3(),saveat=0.01,abstol = 1e-6, reltol = 1e-3)
    #sol = solve(prob,Rodas5P(),saveat=1.0,abstol = 1e-6, reltol = 1e-3)
    sol = solve(prob,Rosenbrock23(),saveat=0.1,abstol = 1e-6, reltol = 1e-3)
 
    #= p1=plot!(sol,marker=(:circle))
  savefig(p1, "plot_144_bs3_1e-3_ft10.png") =#
    
end
@time test()
@time odeDiffEquPackage()   =#