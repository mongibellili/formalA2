using formalA2
using XLSX
#using BenchmarkTools
println("after loading")
#include("./typeB1_1.jl")
include("./typeB1.jl")
include("./typeB2.jl")



println("after INCLUDING")

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
    #save_Sol(solnmliqss2)
    solnmliqss2Interp=solInterpolated(solnmliqss2,0.01)
    er1=getError(solnmliqss2Interp,1,x1)  
    er2=getError(solnmliqss2Interp,2,x2) 
    # timenmliqss=@belapsed QSS_Solve($odeprob,nmliqss2(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=100.0)
    resnmliqss= ("$(solnmliqss2.sysName)",(er1+er2)/2,solnmliqss2.totalSteps,solnmliqss2.simulStepCount)
   # @show resnmliqss
end


function mainTestB1(solver::formalA2.QSSAlgorithm{solType, V})where {V,solType} 
    absTol=1e-6
    relTol=1e-3
    ft=2.5   
    cResults=[]
    #for fun in funs
    total=3040
    for k=5261:5260+total
      fn=Symbol("B1_$k")
      fun=getfield(Main, fn)
      if 6000<k<6002 @show k end
      if 7000<k<7002 @show k end
      if 8000<k<8002 @show k end
     
    #  @show fun
      push!(cResults,solveProblem(fun,ft,solver,absTol,relTol))
    end
    #@show cResults
    XLSX.openxlsx("LTI_B1_2_$(solType)_$(V).xlsx", mode="w") do xf
      sheet = xf[1]
      sheet["A1"] = "LTI_B"
      sheet["A2"] = "$(solType)_relTol=$relTol"
      sheet["A3"] = collect(("problem","error","totalSteps","simul_steps"))
      for i=4:total+3
        sheet["A$i"] = collect(cResults[i-3])
      end
    end
end

mainTestB1(nmliqss1())#normal analytic

function mainTestB2(solver::formalA2.QSSAlgorithm{solType, V})where {V,solType} 
  absTol=1e-6
  relTol=1e-3
  ft=2.5   
  cResults=[]
  #for fun in funs
  total=5160
  for k=10361:total+10360
    fn=Symbol("B2_$k")
    fun=getfield(Main, fn)
  #  @show fun
    push!(cResults,solveProblem(fun,ft,solver,absTol,relTol))
  end
  #@show cResults
  XLSX.openxlsx("LTI_B2_$(solType)_$(V).xlsx", mode="w") do xf
    sheet = xf[1]
    sheet["A1"] = "LTI_B"
    sheet["A2"] = "$(solType)_relTol=$relTol"
    sheet["A3"] = collect(("problem","error","totalSteps","simul_steps"))
    for i=4:total+3
      sheet["A$i"] = collect(cResults[i-3])
    end
  end
end

mainTestB2(nmliqss1())#normal analytic