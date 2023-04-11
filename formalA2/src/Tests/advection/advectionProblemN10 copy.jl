using formalA2
#using formalqss
using XLSX
using BSON
#using BenchmarkTools
#using Plots

#using OrdinaryDiffEq
include("/home/unknown/models/Advection.jl")
using BSON
#include("D://Advection.jl") 
function test()
    
    BSON.@load "bson_base/solVectAdvection_N10_Feagin14e-12.bson" solFeagin14VectorN10
    odeprob = @NLodeProblemLoop begin
        u =[1.0, 1.0, 1.0, 0.0, 0.0,0.0,0.0,0.0,0.0,0.0]
        discrete = [0.0]
        du[1] = -120.0*u[1]+10.0*u[2]+1000.0*u[1]*u[1]*(1.0-u[1]) 
         
        for i in 2:9  
            du[i]=10.0*u[i+1]+110.0*u[i-1]-120.0*u[i]+1000.0*u[i]*u[i]*(1.0-u[i]) ;
       end 

        du[10]=120.0*u[9]-120.0*u[10]+1000.0*u[10]*u[10]*(1.0-u[10]) 
        
    end

  #=   solmliqss2=QSS_Solve_from_model(N10,odeprob,10.0,liqss2(),saveat(0.01),0.0,1e-5,1e-5)
    solmliqss2Interp=solInterpolated(solmliqss2,0.01,10.0)
    err=getAverageErrorByRodas(solFeagin14VectorN10,solmliqss2Interp)
    liqssRes= ("liqss",err,solmliqss2.totalSteps,solmliqss2.simulStepCount)
    #save_Sol(solmliqss2,"x1",1;xlims=(0.0,10.0),ylims=(0.87585,0.8759))
   #=  save_Sol(solmliqss2,"x2",2;xlims=(0.0,10.0),ylims=(0.98752,0.9876))
    save_Sol(solmliqss2,"x3 ",3;xlims=(0.0,10.0),ylims=(0.99875,0.9988))
    save_Sol(solmliqss2,"x4 ",4;xlims=(0.0,10.0),ylims=(0.99987,0.99990))
    save_Sol(solmliqss2,"x5-x10 ",5,6,7,8,9,10;xlims=(0.0,10.0),ylims=(0.9999855,1.0001))
 =#
    solmliqss2=QSS_Solve_from_model(N10,odeprob,10.0,mliqss2(),saveat(0.01),0.0,1e-5,1e-5)
    solmliqss2Interp=solInterpolated(solmliqss2,0.01,10.0)
    err=getAverageErrorByRodas(solFeagin14VectorN10,solmliqss2Interp)
    mliqssRes= ("mliqss",err,solmliqss2.totalSteps,solmliqss2.simulStepCount)
   # save_Sol(solmliqss2,"x1",1;xlims=(0.0,10.0),ylims=(0.87585,0.8759))
    #= save_Sol(solmliqss2,"x2",2;xlims=(0.0,10.0),ylims=(0.98752,0.9876))
    save_Sol(solmliqss2,"x3 ",3;xlims=(0.0,10.0),ylims=(0.99875,0.9988))
    save_Sol(solmliqss2,"x4 ",4;xlims=(0.0,10.0),ylims=(0.99987,0.99990))
    save_Sol(solmliqss2,"x5-x6 ",5,6;xlims=(0.0,10.0),ylims=(0.99997,1.00001))
    save_Sol(solmliqss2,"x7-x10 ",7,8,9,10;xlims=(0.0,10.0),ylims=(0.9999855,1.0001)) =#

    solmliqss2=QSS_Solve_from_model(N10,odeprob,10.0,nliqss2(),saveat(0.01),0.0,1e-5,1e-5)
    solmliqss2Interp=solInterpolated(solmliqss2,0.01,10.0)
    err=getAverageErrorByRodas(solFeagin14VectorN10,solmliqss2Interp)
    nliqssRes= ("nliqss",err,solmliqss2.totalSteps,solmliqss2.simulStepCount) =#
   # save_Sol(solmliqss2,"x1",1;xlims=(0.0,10.0),ylims=(0.87585,0.8759))
#=     save_Sol(solmliqss2,"x2",2;xlims=(0.0,10.0),ylims=(0.98752,0.9876))
    save_Sol(solmliqss2,"x3 ",3;xlims=(0.0,10.0),ylims=(0.99875,0.9988))
    save_Sol(solmliqss2,"x4 ",4;xlims=(0.0,10.0),ylims=(0.99987,0.99990))
    save_Sol(solmliqss2,"x5-x10 ",5,6,7,8,9,10;xlims=(0.0,10.0),ylims=(0.9999855,1.0001)) =#

    solmliqss2=QSS_Solve_from_model(N10,odeprob,10.0,nmliqss2(),saveat(0.01),0.0,1e-5,1e-5)
    solmliqss2Interp=solInterpolated(solmliqss2,0.01,10.0)
    err=getAverageErrorByRodas(solFeagin14VectorN10,solmliqss2Interp)
    nmliqssRes= ("nmliqss",err,solmliqss2.totalSteps,solmliqss2.simulStepCount)
    @show nmliqssRes
   # save_Sol(solmliqss2,"x1",1;xlims=(0.0,10.0),ylims=(0.87585,0.8759))
#=     save_Sol(solmliqss2,"x2",2;xlims=(0.0,10.0),ylims=(0.98752,0.9876))
    save_Sol(solmliqss2,"x3 ",3;xlims=(0.0,10.0),ylims=(0.99875,0.9988))
    save_Sol(solmliqss2,"x4 ",4;xlims=(0.0,10.0),ylims=(0.99987,0.99990))
    save_Sol(solmliqss2,"x5-x10 ",5,6,7,8,9,10;xlims=(0.0,10.0),ylims=(0.9999855,1.0001)) =#
  #=   XLSX.openxlsx("sysN10 all solvers___uii=x-ujj=x-___.xlsx", mode="w") do xf
        sheet = xf[1]
        sheet["A1"] = "sysN10 all solvers_"
        sheet["A4"] = collect(("solver","error","totalSteps","simul_steps"))
        sheet["A5"] = collect(liqssRes)
        sheet["A6"] = collect(mliqssRes)
        sheet["A7"] = collect(nliqssRes)
        sheet["A8"] = collect(nmliqssRes)
   end   =#
end
#@btime 
test()


