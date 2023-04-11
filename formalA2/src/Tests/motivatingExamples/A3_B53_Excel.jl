
using XLSX
using formalA2
#using formalqss
#using formalqssAB
using BenchmarkTools
#= using StaticArrays
using BSON =#
include("/home/unknown/models/BProblem.jl")
include("/home/unknown/models/classicProblem.jl")

function test3(ft::Float64,solver::Val{V})where {V} 
     odeprob = @NLodeProblem begin
           #------------twoVarSys13----------
           u = [-4.0, 4.0]
           discrete = [0.0]
           du[1] = -u[1]-u[2]+0.2
           du[2] =u[1]-u[2]+1.2    
     end
     solmliqss2=QSS_Solve_from_model(twoVarSys13,odeprob,ft,solver,saveat(0.01),0.0,1e-6,1e-3)
     solmliqss2Interp=solInterpolated(solmliqss2,0.1,100.0)
    # save_Sol(solmliqss2,"x1",1;xlims=(30.0,40.0),ylims=(-0.55,-0.45))
     u1, u2 = 0.0 - 1.0im, -0.0 + 1.0im
    λ1, λ2 = -1.0 - 1.0im, -1.0 + 1.0im
   
    c1, c2 = 1.65 - 1.75im, 1.65 + 1.75im
    xp1, xp2 = -0.5, 0.7
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
   # saveAbsoluteError(solmliqss2Interp,1,x1)  
   # saveAbsoluteError(solmliqss2Interp,2,x2) 
    er1=Float64(getError(solmliqss2Interp,1,x1))  
    er2=Float64(getError(solmliqss2Interp,2,x2)) 
    ("$(solmliqss2.algName)",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
end

function test53(ft::Float64,solver::Val{V})where {V} 
     # BSON.@load "sysB/relaxedAE_sysB53_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
     solmliqss2=QSS_Solve_from_model(B53,probB53(),ft,solver,saveat(0.01),0.0,1e-6,1e-3)
     solmliqss2Interp=solInterpolated(solmliqss2,0.1,100.0)
     #save_Sol(solmliqss2,"x1",1;xlims=(30.0,40.0),ylims=(-0.000005,0.000005))
     u1, u2 = -8.73522174738572, -7.385745994549763
     λ1, λ2 = -10.841674966758294, -9.168325033241706
     
     c1, c2 = 121.14809142478035, -143.14809142478035
     xp1, xp2 = 0.0, 20.0
     x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
     x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
     
     er1=getError(solmliqss2Interp,1,x1)  
     er2=getError(solmliqss2Interp,2,x2) 
          ("$(solmliqss2.algName)",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
  end


  function test13(ft::Float64,solver::Val{V})where {V} 
     #BSON.@load "sysB/relaxedAE_sysB13_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
     solmliqss2=QSS_Solve_from_model(B13,probB13(),ft,solver,saveat(0.01),0.0,1e-6,1e-3)
     solmliqss2Interp=solInterpolated(solmliqss2,0.1,100.0)
     #save_Sol(solmliqss2,"x1",1;xlims=(30.0,35.0),ylims=(-0.0005,0.0005))
     u1, u2 = -0.06805253887502617, -0.5248046039821167
     λ1, λ2 = -3.523632227874817, -0.32636777212518275
    
     c1, c2 = -27.467202112518294, 5.467202112518295
     xp1, xp2 = 3.552713678800501e-15, 20.0
     x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
     x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
  
  
     er1=getError(solmliqss2Interp,1,x1)  
     er2=getError(solmliqss2Interp,2,x2) 
    
      ("$(solmliqss2.algName)",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
    
 end


function globaltest3()
          function liqsstestOneSys()
               ft=100.0
               solver=liqss2()
               test3(ft,solver)
          end
          #@time
          resultliqss2=liqsstestOneSys()
          function mliqsstestOneSys()
               ft=100.0
               solver=mliqss2()
               test3(ft,solver)
          end
          resultmliqss2=mliqsstestOneSys()
          #@show resultmliqss2
          function nliqsstestOneSys()
               ft=100.0
               solver=nliqss2()
               test3(ft,solver) 
          end
          #@time
          resultnliqss2=nliqsstestOneSys()
          function nmliqsstestOneSys()
               ft=100.0
               solver=nmliqss2()
               test3(ft,solver)
          end
          resultnmliqss2=nmliqsstestOneSys()
          XLSX.openxlsx("sys3 all solvers______.xlsx", mode="w") do xf
               sheet = xf[1]
               sheet["A1"] = "sys3 all solvers_"
               sheet["A4"] = collect(("solver","error","totalSteps","simul_steps"))
               sheet["A5"] = collect(resultliqss2)
               sheet["A6"] = collect(resultmliqss2)
               sheet["A7"] = collect(resultnliqss2)
               sheet["A8"] = collect(resultnmliqss2)
          end  
end
globaltest3()
function globaltestB13()
     function liqsstestOneSys()
          ft=100.0
          solver=liqss2()
          test13(ft,solver)
     end 
     #@time
     resultliqss2=liqsstestOneSys()
     function mliqsstestOneSys()
          ft=100.0
          solver=mliqss2()
          test13(ft,solver)
     end
     resultmliqss2=mliqsstestOneSys()
     #@show resultmliqss2
     function nliqsstestOneSys()
          ft=100.0
          solver=nliqss2()
          test13(ft,solver) 
     end
     #@time
     resultnliqss2=nliqsstestOneSys()
     function nmliqsstestOneSys()
          ft=100.0
          solver=nmliqss2()
          test13(ft,solver)
     end
     resultnmliqss2=nmliqsstestOneSys()
     #@show resultnmliqss2
 

     XLSX.openxlsx("sys13 all solvers______.xlsx", mode="w") do xf
          sheet = xf[1]
          sheet["A1"] = "sys13 all solvers_"
          sheet["A4"] = collect(("solver","error","totalSteps","simul_steps"))
          sheet["A5"] = collect(resultliqss2)
          sheet["A6"] = collect(resultmliqss2)
          sheet["A7"] = collect(resultnliqss2)
          sheet["A8"] = collect(resultnmliqss2)
          end  

end
globaltestB13()
function globaltestB53()
     function liqsstestOneSys()
          ft=100.0
          solver=liqss2()
          test53(ft,solver)
     end 
     #@time
     resultliqss2=liqsstestOneSys()
     function mliqsstestOneSys()
          ft=100.0
          solver=mliqss2()
          test53(ft,solver)
     end
     resultmliqss2=mliqsstestOneSys()
     #@show resultmliqss2
     function nliqsstestOneSys()
          ft=100.0
          solver=nliqss2()
          test53(ft,solver) 
     end
     #@time
     resultnliqss2=nliqsstestOneSys()
     function nmliqsstestOneSys()
          ft=100.0
          solver=nmliqss2()
          test53(ft,solver)
     end
     resultnmliqss2=nmliqsstestOneSys()
     #@show resultnmliqss2
 

     XLSX.openxlsx("sys53 all solvers______.xlsx", mode="w") do xf
          sheet = xf[1]
          sheet["A1"] = "sys53 all solvers_"
          sheet["A4"] = collect(("solver","error","totalSteps","simul_steps"))
          sheet["A5"] = collect(resultliqss2)
          sheet["A6"] = collect(resultmliqss2)
          sheet["A7"] = collect(resultnliqss2)
          sheet["A8"] = collect(resultnmliqss2)
          end  

end
globaltestB53()