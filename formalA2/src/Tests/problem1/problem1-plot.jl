using formalA2
#using OrdinaryDiffEq
#using BenchmarkTools
#using Plots


#include("D://QS_Solver//qss//src//models//classicProblem.jl") 
#include("/home/unknown/formalqssAB/formalqssAB/src/models/classicProblem.jl")# where you saved the model
#include("/home/unknown/formalqssAB/formalqssAB/src/models/InvestigateBProblem.jl") 
include("/home/unknown/formalA2/formalA2/src/models/NonLinear.jl") 

function test()
    odeprob = @NLodeProblem begin
      ##########lorenz equation
      #= u = [1.0, 0.0,0.0]
      discrete = [0.0]
      du[1] = 10.0*(u[2]-u[1])
      du[2] = u[1]*(28.0-u[3]) - u[2]
      du[3] = u[1]*u[2] - (2.66667)*u[3] =#
      #######the Lotka-Volterra model...
      u = [1.0, 1.0]
      discrete = [0.0]
      du[1] = 2.0 * u[1] - 1.2 * u[1]*u[2]
      du[2] = -3.0 * u[2] + u[1]*u[2]

      ###########Van der pol
      #= u = [0.0, 1.7320]
      discrete = [0.0]
      du[1] = u[2]
      du[2] = (1.0-u[1]*u[1]) * u[2] - u[1] =#
      #############oregonator
     #=  u = [1.0, 0.0,0.0]
      discrete = [0.0]
      du[1] = 100.8*(   9.523809523809524e-5 * u[2] - u[1] * u[2] + u[1] * ( 1.0 - u[1] ) ) 
      du[2] = 40320.0*( - 9.523809523809524e-5 * u[2] - u[1] * u[2] +  u[3] ) 
      du[3] = u[1] - u[3] =#

  end
  solqss2=QSS_Solve_from_model(loktavoltera,odeprob,10.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
  save_Sol(solqss2,"formalAE")
  solqss2=QSS_Solve_from_model(loktavoltera,odeprob,10.0,mliqss3(),saveat(0.1),0.0,1e-6,1e-3)
  save_Sol(solqss2,"formalAE")
  

end

#@time 
test()

function test2()
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
    u = [0.0, 1.7320]
    discrete = [0.0]
    du[1] = u[2]
    du[2] = (1.0-u[1]*u[1]) * u[2] - u[1]
    #############oregonator
   #=  u = [1.0, 0.0,0.0]
    discrete = [0.0]
    du[1] = 100.8*(   9.523809523809524e-5 * u[2] - u[1] * u[2] + u[1] * ( 1.0 - u[1] ) ) 
    du[2] = 40320.0*( - 9.523809523809524e-5 * u[2] - u[1] * u[2] +  u[3] ) 
    du[3] = u[1] - u[3] =#

end
 solqss2=QSS_Solve_from_model(vanderpol,odeprob,10.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
save_Sol(solqss2,"formalAE")
solqss2=QSS_Solve_from_model(vanderpol,odeprob,10.0,mliqss3(),saveat(0.1),0.0,1e-6,1e-3)
save_Sol(solqss2,"formalAE")



end

#@time 
test2()