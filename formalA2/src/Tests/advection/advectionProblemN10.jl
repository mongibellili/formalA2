using formalA2
#using formalqssA
#using BenchmarkTools
#using Plots
#using OrdinaryDiffEq

include("/home/unknown/models/Advection.jl")
using BSON
#include("D://Advection.jl") 
function test()
    odeprob = @NLodeProblemLoop begin
        u =[1.0, 1.0, 1.0, 0.0, 0.0,0.0,0.0,0.0,0.0,0.0]
        discrete = [0.0]
        du[1] = -120.0*u[1]+10.0*u[2]+1000.0*u[1]*u[1]*(1.0-u[1]) 
         
        for i in 2:9  
            du[i]=10.0*u[i+1]+110.0*u[i-1]-120.0*u[i]+1000.0*u[i]*u[i]*(1.0-u[i]) ;
       end 

        du[10]=120.0*u[9]-120.0*u[10]+1000.0*u[10]*u[10]*(1.0-u[10]) 
        
    end

    # solmliqss2=QSS_Solve(odeprob,10.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
    #save_prob_to_model(odeprob,"d://Advection.jl","N10") #any location you want   

    solmliqss2=QSS_Solve_from_model(N10,odeprob,10.0,mliqss2(),saveat(0.01),0.0,1e-3,1e-3)
    #save_Sol(solmliqss2,"formalAEM",2,3,4,10;xlims=(0.0,1.0),ylims=(0.998,1.002))
    #= solmliqss2Interp=solInterpolated(solmliqss2,0.01,0.6)
    BSON.@save "solAdvection_N10_mliqss2e-6ft06.bson" solmliqss2Interp =#
    save_Sol(solmliqss2,"formalAEM_flags",3,5,6,10;xlims=(0.0,10.0),ylims=(0.995,1.005))
  # save_SolVar(solmliqss2,1000)
  
 #save_Sol(solmliqss2)
end
#@btime 
test()


