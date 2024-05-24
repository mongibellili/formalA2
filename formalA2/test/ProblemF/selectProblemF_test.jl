using formalA2
using XLSX
#using BenchmarkTools
include("typeF.jl")
function solveProblem(prblem::Function,ft::Float64,solver::formalA2.QSSAlgorithm{solType, V},absTol,relTol)where {V,solType} 
    pr=prblem()
    prob=pr[1]
    x1=pr[2]
    x2=pr[3]
    timenmliqss=0.0
    #= absTol=1e-5
    relTol=1e-3 =#
    tspan=(0.0,ft)
    solnmliqss2=solve(prob,solver,abstol=absTol,saveat=0.01,reltol=relTol,tspan)
   # save_Sol(solnmliqss2)
    solnmliqss2Interp=solInterpolated(solnmliqss2,0.01)
    er1=getError(solnmliqss2Interp,1,x1)  
    er2=getError(solnmliqss2Interp,2,x2) 
    # timenmliqss=@belapsed QSS_Solve($odeprob,nmliqss2(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=100.0)
    resnmliqss= ("$(solnmliqss2.sysName)",(er1+er2)/2,solnmliqss2.totalSteps,solnmliqss2.simulStepCount,timenmliqss)
   # @show resnmliqss
end


function mainTest(solver::formalA2.QSSAlgorithm{solType, V})where {V,solType} 
    absTol=1e-6
    relTol=1e-3
    ft=100.0
    
    #funs=[F_1,F_2,F_3,F_4,F_5,F_6,F_7,F_8,F_9,F_10,F_11,F_12,F_13,F_14,F_15,F_16,F_17,F_18,F_19,F_20,F_21,F_22,F_23,F_24,F_25,F_26,F_27,F_28,F_29,F_30,F_31,F_32,F_33,F_34,F_35,F_36,F_37,F_38,F_39,F_40,F_41,F_42,F_43,F_44,F_45,F_46,F_47,F_48,F_49,F_50,F_51,F_52,F_53,F_54,F_55,F_56,F_57,F_58,F_59,F_60,]
    funs=[F_61]
    results=[]
    for fun in funs
      push!(results,solveProblem(fun,ft,solver,absTol,relTol))
    end
    #@show results
    XLSX.openxlsx("LTI_F_$(solType).xlsx", mode="w") do xf
      sheet = xf[1]
      sheet["A1"] = "LTI_F"
      sheet["A2"] = "$(solType)_relTol=$relTol"
      sheet["A3"] = collect(("problem","error","totalSteps","simul_steps","time"))
      for i=4:length(funs)+3
        sheet["A$i"] = collect(results[i-3])
      end
    end
end


#mainTest(nmliqss1())#normal analytic
#mainTest(mliqss1())#golden search analytic
#mainTest(nliqss1()) #iters
mainTest(mliqssBounds1()) #analyticBounds