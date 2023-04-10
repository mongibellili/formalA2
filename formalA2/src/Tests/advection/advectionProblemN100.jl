using formalA2
#using formalqssAC
#using formalqssA
#using BenchmarkTools
using BSON
#using Plots
#using OrdinaryDiffEq

#include("D://QS_Solver//qss//src//models//classicProblem.jl") 
#D:\qssv01\qssv0.1\qssv01\src\models\classicProblem.jl
#include("D://Advection.jl") 
include("/home/unknown/models/Advection.jl")
function test()
    odeprob = @NLodeProblemLoop begin
      u = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        discrete = [0.0]
        du[1] = -120.0*u[1]+10.0*u[2]+1000.0*u[1]*u[1]*(1.0-u[1]) 
         
        for i in 2:99  
            du[i]=10.0*u[i+1]+110.0*u[i-1]-120.0*u[i]+1000.0*u[i]*u[i]*(1.0-u[i]) ;
       end 

        du[100]=120.0*u[99]-120.0*u[100]+1000.0*u[100]*u[100]*(1.0-u[100]) 
        
    end

    # solmliqss2=QSS_Solve(odeprob,10.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
   # save_prob_to_model(odeprob,"/home/unknown/models/Advection.jl","N100") 

    solmliqss2=QSS_Solve_from_model(N100,odeprob,10.0,mliqss2(),saveat(0.01),0.0,1e-5,1e-5)
    solmliqss2Interp=solInterpolated(solmliqss2,0.01,10.0)
    BSON.@save "solAdvection_N100_NewLiqss-5.bson" solmliqss2Interp
    @show solmliqss2Interp.totalSteps

  # save_SolVar(solmliqss2,1000)
  save_Sol(solmliqss2,"NewLiqss2_N100",3,4,5,6;xlims=(0.0,10.0),ylims=(0.99995,1.00002))
 #save_Sol(solmliqss2)
end
#@time 
test()
