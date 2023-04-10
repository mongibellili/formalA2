using formalA2
#using formalqssAC
#using formalqssA
using BenchmarkTools
#using Plots
#using OrdinaryDiffEq

include("/home/unknown/models/Advection.jl")
function test()
    odeprob = @NLodeProblemLoop begin
        u = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        discrete = [0.0]
        du[1] = -120.0*u[1]+10.0*u[2]+1000.0*u[1]*u[1]*(1.0-u[1]) 
         
        for i in 2:49  
            du[i]=10.0*u[i+1]+120.0*u[i-1]-120.0*u[i]+1000.0*u[i]*u[i]*(1.0-u[i]) ;
       end 
        du[50]=120.0*u[49]-120.0*u[50]+1000.0*u[50]*u[50]*(1.0-u[50]) 
        
    end

   # solmliqss2=QSS_Solve(odeprob,10.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
   solmliqss2=QSS_Solve_from_model(N50,odeprob,10.0,mliqss3(),saveat(0.1),0.0,1e-6,1e-3)
   # solliqss2=QSS_Solve_from_model(twoVarSys14,odeprob,10.0,qss3(),saveat(0.1),0.0,1e-6,1e-3)
  # save_SolVar(solmliqss2,1000)
   #save_SolVar(solmliqss2,400,600,800,1000)
  # save_Sol(solmliqss2,3,10,5,6)
 #  save_SolVar(solmliqss2,1000)
 save_Sol(solmliqss2)
end
#@time 
test()
