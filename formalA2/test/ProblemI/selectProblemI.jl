using formalA2
using XLSX
using Plots
#using BenchmarkTools

#include("C:/Users/belli/rel/formalA2/formalA2/test/typeH.jl")
include("./typeI.jl")
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
    solnmliqss2Interp=solInterpolated(solnmliqss2,0.01)
   #=  p=getPlot!(solnmliqss2)
    p=plot!(x1)
    p=plot!(x2)
    savefig(p, "plotInterpl_$(solType)_$V.png") =#
   
    er1=getError(solnmliqss2Interp,1,x1)  
    er2=getError(solnmliqss2Interp,2,x2) 
    # timenmliqss=@belapsed QSS_Solve($odeprob,nmliqss2(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=100.0)
    resnmliqss= ("$(solnmliqss2.sysName)",(er1+er2)/2,solnmliqss2.totalSteps,solnmliqss2.simulStepCount,timenmliqss)
   # @show resnmliqss
end


function mainTest(solver::formalA2.QSSAlgorithm{solType, V})where {V,solType} 
    absTol=1e-5
    relTol=1e-3
    ft=3.0
    
  funs=[I_37,I_38,I_39,I_40,I_41,I_42,I_43,I_44,I_45,I_46,I_47,I_48,I_49,I_50,I_51,I_52,I_53,I_54,I_55,I_56,I_57,I_58,I_59,I_60,I_61,I_62,I_63,I_64,I_65,I_66,I_67,I_68,I_69,I_70,I_71,I_72,I_73,I_74,I_75,I_76,I_77,I_78,I_79,I_80,I_81,I_82,I_83,I_84,I_85,I_86,I_87,I_88,I_89,I_90,I_91,I_92,I_93,I_94,I_95,I_96,I_97,I_98,I_99,I_100,I_101,I_102,I_103,I_104,I_105,I_106,I_107,I_108,I_109,I_110,I_111,I_112,I_113,I_114,I_115,I_116,I_117,I_118,I_119,I_120,I_121,I_122,I_123,I_124,I_125,I_126,I_127,I_128,I_129,I_130,I_131,I_132,I_133,I_134,I_135,I_136,I_137,I_138,I_139,I_140,I_141,I_142,I_143,I_144,I_145,I_146,I_147,I_148,I_149,I_150,I_151,I_152,I_153,I_154,I_155,I_156,I_157,I_158,I_159,I_160,I_161,I_162,I_163,I_164,I_165,I_166,I_167,I_168,I_169,I_170,I_171,I_172,I_173,I_174,I_175,I_176,I_177,I_178,I_179,I_180,]

    results=[]
    cntr=0
    for fun in funs
      cntr+=1
      if 100<cntr <102 @show fun end
      if 200<cntr <202 @show fun end
      if 300<cntr <302 @show fun end
      if 400<cntr <402 @show fun end
      push!(results,solveProblem(fun,ft,solver,absTol,relTol))
    end
    
   # @show results
    XLSX.openxlsx("LTI_I_$(solType).xlsx", mode="w") do xf
      sheet = xf[1]
      sheet["A1"] = "LTI_I"
      sheet["A2"] = "$(solType)_relTol=$relTol"
      sheet["A3"] = collect(("problem","error","totalSteps","simul_steps","time"))
      for i=4:length(funs)+3
        sheet["A$i"] = collect(results[i-3])
      end
    end
end


mainTest(nmliqss1())#normal analytic
#mainTest(nmliqss2())#normal analytic
#mainTest(mliqss1())#golden search analytic
#mainTest(nliqss1()) #iters
#mainTest(mliqssBounds1()) #analyticBounds