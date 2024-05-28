
using formalA2
#using XLSX
#using BenchmarkTools
using BSON
#using TimerOutputs
#using Plots
function test(case,solvr)
   #  BSON.@load "formalA2/ref_bson/solVectlorenz_Feagin14_e-8.bson"  sollorenzFeagin14
  absTol=1e-6
     relTol=1e-3
      odeprob = @NLodeProblem begin
         #sys b53
         name=(lorenz,)
         u = [1.0, 1.0,0.0]
         du[1] = 10.0*(u[2]-u[1])
         du[2] =u[1]*(28.0-u[3])-u[2]
         du[3] =u[2]*u[1]-2.66667*u[3]
     end  


     timenmliqss=0.0;er1=0.0;er2=0.0
     
     println("start LTI solving")
     tspan=(0.0,10.0)
     solnmliqss=solve(odeprob,solvr,abstol=absTol,saveat=0.01,reltol=relTol,tspan#= ,maxErr=10000*relTol =#)
     save_Sol(solnmliqss) 
    # @show solnmliqss.totalSteps
   # save_Sol(solnmliqss,note="xi10dxi")

  # solnmliqssInterp=solInterpolated(solnmliqss,0.01)
   #@show solnmliqssInterp.savedVars
   err=0.0
 # err=getAverageErrorByRefs( sollorenzFeagin14,solnmliqssInterp)
#  @show err4,solnmliqss.totalSteps
   # timenmliqss=@belapsed solve($odeprob,$solvr,abstol=$absTol,saveat=0.01,reltol=$relTol,tspan#= ,maxErr=1000*$relTol =#)
     resnmliqss1E_2= ("$(solnmliqss.algName)",relTol,err,solnmliqss.totalSteps,solnmliqss.simulStepCount,timenmliqss)
     @show resnmliqss1E_2

    
 
end

case="order1_"
test(case,nmliqss1())
#test(case,mliqssBounds1())
#test(case,nliqss1())