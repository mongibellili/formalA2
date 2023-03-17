using relaxedqssA
#= using RuntimeGeneratedFunctions
using StaticArrays
using SymEngine
using ExprTools
using MacroTools: postwalk, prewalk, @capture
#= import Base.:-
import Base.:+
import Base.:* =#
using Plots: plot!,plot,savefig
using Dates: now,year,month,day,hour,minute,second #fortimestamp =#


using BSON
include("/home/unknown/formalqssA/formalqssA/src/models/classicProblem.jl") 
function test11()
    odeprob = @NLodeProblem begin
         #sys b11
         u = [1.0, 0.0]
       discrete = [0.0]
       du[1] = 0.01*u[2]
       du[2] =-100.0*u[1]-100.0*u[2]+2020.0

     
    end
  display(odeprob)
 # BSON.@save "sys12.bson" odeprob

 #=  BSON.@load  "sys12.bson" odeprob
   solmliqss2=QSS_Solve_from_model(twoVarSys12,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
   save_Sol(solmliqss2) =#
   #=  save_Sol(solmliqss2)
    BSON.@save "sys1/relaxedqssA_sys12_mliqss2_errorE-3_step0.01_ft1000.bson" solmliqss2
    solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
    BSON.@save "sys1/relaxedqssA_sys12_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp =#
   
end
test11() 
