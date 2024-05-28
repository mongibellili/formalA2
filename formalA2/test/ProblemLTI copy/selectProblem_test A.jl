using formalA2
using XLSX
using Plots
#using BenchmarkTools

#= include("./typeA1_1.jl")
include("./typeAA1.jl")
include("./typeA1_2.jl")
include("./typeA2_1.jl")
include("./typeA2_2.jl")
include("./typeA2_3.jl")
include("./typeAA2_1.jl")
include("./typeAA2_2.jl") =#

function A2_839() 
  odeprob = @NLodeProblem begin
                 name=(A2_839,)
                 u = [-1.0, -2.0]
                 du[1] = -0.75*u[1]+-21.0*u[2]+-20.0
                 du[2] =-2.0*u[1]+-21.0*u[2]+-20.0
             end  
  x1(t)=-1.0481691867790153*0.9482325884620754*exp(-22.89646517692415*t)+0.0005501391599678998*-11.073232588462075*exp(1.1464651769241527*t)+0.0
  x2(t)=-1.0481691867790153*exp(-22.89646517692415*t)+0.0005501391599678998*exp(1.1464651769241527*t)+-0.9523809523809524
  return (odeprob,x1,x2) 
  end 


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
   save_Sol(solnmliqss2)
    solnmliqss2Interp=solInterpolated(solnmliqss2,0.01)
    er1=getError(solnmliqss2Interp,1,x1)  
    er2=getError(solnmliqss2Interp,2,x2) 
   
     @show (er1+er2)/2
    # timenmliqss=@belapsed QSS_Solve($odeprob,nmliqss2(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=100.0)
    resnmliqss= ("$(solnmliqss2.sysName)",solnmliqss2.totalSteps,solnmliqss2.simulStepCount)
   # @show resnmliqss
end


#= function mainTestA1(solver::formalA2.QSSAlgorithm{solType, V})where {V,solType} 
  absTol=1e-6
  relTol=1e-3
  ft=100.0   
  cResults=[]
  #for fun in funs
  total=5800
  for k=1:total
    fn=Symbol("A1_$k")
    fun=getfield(Main, fn)
    if 500<k<502 @show k end
    if 1000<k<1002 @show k end
    if 2000<k<2002 @show k end
    if 3000<k<3002 @show k end
    if 4000<k<4002 @show k end
    
  #  @show fun
    push!(cResults,solveProblem(fun,ft,solver,absTol,relTol))
  end
  #@show cResults
  XLSX.openxlsx("LTI_ft100_A1_$(solType)_$(V).xlsx", mode="w") do xf
    sheet = xf[1]
    sheet["A1"] = "LTI_A1"
    sheet["A2"] = "$(solType)_relTol=$relTol"
    sheet["A3"] = collect(("problem","totalSteps","simul_steps"))
    for i=4:total+3
      sheet["A$i"] = collect(cResults[i-3])
    end
  end
end

mainTestA1(nmliqss1())#normal analytic
 =#
function mainTestA2(solver::formalA2.QSSAlgorithm{solType, V})where {V,solType} 
  absTol=1e-6
  relTol=1e-3
  ft=2.5  
  cResults=[]
  #for fun in funs
  total=12540
  for k=839:839
    fn=Symbol("A2_$k")
    fun=getfield(Main, fn)
  #  @show fun
 #=  if 500<k<502 @show k end
  if 1000<k<1002 @show k end
  if 2000<k<2002 @show k end
  if 3000<k<3002 @show k end
  if 4000<k<4002 @show k end
  if 5000<k<5002 @show k end
  if 6000<k<6002 @show k end
  if 7000<k<7002 @show k end
  if 8000<k<8002 @show k end
  if 9000<k<9002 @show k end
  if 10000<k<10002 @show k end
  if 11000<k<11002 @show k end
  if 12000<k<12002 @show k end =#
    push!(cResults,solveProblem(fun,ft,solver,absTol,relTol))
  end
  #@show cResults
 #=  XLSX.openxlsx("LTI_A2_ft100_$(solType)_$(V).xlsx", mode="w") do xf
    sheet = xf[1]
    sheet["A1"] = "LTI_A2"
    sheet["A2"] = "$(solType)_relTol=$relTol"
    sheet["A3"] = collect(("problem","totalSteps","simul_steps"))
    for i=4:total+3
      sheet["A$i"] = collect(cResults[i-3])
    end
  end =#
end

mainTestA2(nmliqss1())#normal analytic