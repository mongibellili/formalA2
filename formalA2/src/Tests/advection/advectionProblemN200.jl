using formalA2
#using formalqssA
#using BenchmarkTools
#using Plots
#using OrdinaryDiffEq

#include("D://QS_Solver//qss//src//models//classicProblem.jl") 
#D:\qssv01\qssv0.1\qssv01\src\models\classicProblem.jl
#include("D://Advection.jl") 
include("/home/unknown/models/Advection.jl")
function test()
    odeprob = @NLodeProblemLoop begin
      u = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        discrete = [0.0]
        du[1] = -120.0*u[1]+10.0*u[2]+1000.0*u[1]*u[1]*(1.0-u[1]) 
         
        for i in 2:199  
            du[i]=10.0*u[i+1]+110.0*u[i-1]-120.0*u[i]+1000.0*u[i]*u[i]*(1.0-u[i]) ;
       end 

        du[200]=120.0*u[199]-120.0*u[200]+1000.0*u[200]*u[200]*(1.0-u[200]) 
        
    end
   # @show odeprob
    # solmliqss2=QSS_Solve(odeprob,10.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
  # save_prob_to_model(odeprob,"/home/unknown/models/Advection.jl","N200") 

   solmliqss2=QSS_Solve_from_model(N200,odeprob,10.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)

  # save_SolVar(solmliqss2,1000)
 # save_Sol(solmliqss2,"formalAEM",3,5,6,10;xlims=(0.0,10.0),ylims=(0.995,1.005))
 save_Sol(solmliqss2)
end
#@btime 
test()
