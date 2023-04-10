using formalA2

#= using BenchmarkTools
#using Plots
#using OrdinaryDiffEq
using BSON
using StaticArrays =#
#using SymEngine
#include("D://Advection.jl") 
include("/home/unknown/models/Advection.jl")
#D:\qssv01\qssv0.1\qssv01\src\models\classicProblem.jl

function test()
    odeprob = @NLodeProblemLoop begin
        u =[1.0, 0.0, 0.0,0.0]
        discrete = [0.0]
        du[1] = -120.0*u[1]+10.0*u[2]+1000.0*u[1]*u[1]*(1.0-u[1]) 
         
        for i in 2:3  
            du[i]=10.0*u[i+1]+120.0*u[i-1]-120.0*u[i]+1000.0*u[i]*u[i]*(1.0-u[i]) ;
       end 

        du[4]=120.0*u[3]-120.0*u[4]+1000.0*u[4]*u[4]*(1.0-u[4]) 
        
    end


 
   # @show odeprob
   # solmliqss2=QSS_Solve(odeprob[1],odeprob[2],10.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3,1e-1)
    #save_prob_to_model(odeprob,"d://Advection.jl","N4") #any location you want   

    solmliqss2=QSS_Solve_from_model(N4,odeprob,10.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)

  # save_SolVar(solmliqss2,1000)
 # save_Sol(solmliqss2,"formalA",3,1;xlims=(0.0,10.0),ylims=(0.0,1.1))
 save_Sol(solmliqss2)
end
#@time
 test()
