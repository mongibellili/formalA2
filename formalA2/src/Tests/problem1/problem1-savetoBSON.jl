using formalqssAB
#using BenchmarkTools
#using Plots
using BSON
#using OrdinaryDiffEq

#include("D://QS_Solver//qss//src//models//classicProblem.jl") 
include("/home/unknown/formalqssAB/formalqssAB/src/models/classicProblem.jl")# where you saved the model
function test()
    odeprob = @NLodeProblem begin
        
        #----------twoVarSys1-------------
        
       #=  u = [1.0, 0.0]
        discrete = [0.0]
        du[1] = u[2]
        du[2] =-u[1]-u[2]  =#
       
        #------------twoVarSys12-----------
        
     #=    u = [1.0, 0.0]
        discrete = [0.0]
        du[1] = 0.01*u[2]
        du[2] =-100.0*u[1]-100.0*u[2]+2020.0 =#

        #------------twoVarSys13----------
        u = [-4.0, 4.0]
        discrete = [0.0]
        du[1] = -u[1]-u[2]+0.2
        du[2] =u[1]-u[2]+1.2
        #------------twoVarSys14----------
        #= u = [-4.0, 4.0]
        discrete = [0.0]
        du[1] = -u[1]-10.0*u[2]+0.2
        du[2] =10.0*u[1]-u[2]+1.2 =#
        #------------twoVarSys15----------
       #=  u = [100.0, 0.0]
        discrete = [0.0]
        du[1] = -4.0*u[1]+0.5*u[2]+1.0
        du[2] =u[1]-0.25*u[2]+0.25 =#
        #------------twoVarSys16----------
       #=  u = [1.0, 2.3]
        discrete = [0.0]
        du[1] = -1.01*u[1]-100.0*u[2]+1000.2
        du[2] =10.1*u[1]-1.001*u[2]+0.2 =#
      
    end
    

    




  # solqss1=QSS_Solve_from_model(twoVarSys16,odeprob,20.0,liqss1(),saveat(0.1),0.0,1e-9,1e-4)    
   #solqss2=QSS_Solve_from_model(twoVarSys13,odeprob,80.0,liqss2(),saveat(0.01),0.0,1e-9,1e-4) 
   # solqss3=QSS_Solve_from_model(twoVarSys1,odeprob,5.0,qss3(),saveat(0.1),0.0,1e-6,1e-3)
   
    #solliqss1=QSS_Solve_from_model(twoVarSys15,odeprob,50.0,liqss1(),saveat(0.5),0.0,1e-9,1e-5)
  #  solliqss2=QSS_Solve_from_model(twoVarSys16,odeprob,100.0,liqss2(),saveat(0.1),0.0,1e-6,1e-3)
   # solliqss3=QSS_Solve_from_model(twoVarSys12,odeprob,100.0,liqss2(),saveat(0.1),0.0,1e-6,1e-3)

    #solmliqss1=QSS_Solve_from_model(twoVarSys15,odeprob,50.0,mliqss1(),saveat(0.5),0.0,1e-9,1e-5)
    solmliqss2=QSS_Solve_from_model(twoVarSys13,odeprob,1000.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
   # solmliqss1=QSS_Solve_from_model(twoVarSys13,odeprob,50.0,mliqss1(),saveat(0.01),0.0,1e-8,1e-4)
    

  #=   u1=(-sqrt(257)-15)/8
    u2=(sqrt(257)-15)/8
    λ1=(-sqrt(257)-17)/8
    λ2=(sqrt(257)-17)/8
    c2=(397-2*(sqrt(257)+15))/sqrt(257)+
    c1=-4-c2
    b1=1.0
    b2=0.25
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+0.5*b1+b2
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+2*b1+8*b2
    plotAbsoluteError(solqss2,1,x1)
    plotAbsoluteError(solmliqss2,1,x1)
    plotRelativeError(solqss2,1,x1)
    plotRelativeError(solmliqss2,1,x1)  =#

    #= plotAbsoluteError(solqss2,2,x2)
    plotAbsoluteError(solmliqss2,2,x2)
    plotRelativeError(solqss2,2,x2)
    plotRelativeError(solmliqss2,2,x2)=#




    # solliqss2=QSS_Solve(odeprob,20.0,qss2())
  # sol=QSS_Solve(odeprob,50.0,qss3(),saveat(0.1),0.0,1e-6,1e-3)
  # display(sol(1,14.0))
 #plotSol(sol)
 #plotSol(solqss1)
 #plotSol(solqss2)
 #plotSol(solliqss1)
 #plotSol(solmliqss1)
 #sol=solInterpolated(solliqss1,1,0.1,1000.0)
 #display(sol)
 #display(plot!(sol[1],sol[2]))
 #sol2=solInterpolated(solliqss1,2,0.1,1000.0)
 BSON.@save "formalAB_sys3_mliqss2_errorE-3_step0.1_ft1000.bson" solmliqss2
 solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
 BSON.@save "formalAB_sys3_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
 #BSON.@save "sys3_mliqss2_errorE-3_step0.1_ft1000_x2.bson" sol2
 #= display(plot!(sol2[1],sol2[2]))
 println("press enter to exit")
 readline()  =#
end
#@time 
test()






#= 
function odeDiffEquPackage()
 
    function funcName(du,u,p,t)# api requires four args
        #= du[1] = -5*u[1]
        du[2] =-u[1]-5*u[2] =#
        #= du[1] = -u[1]-u[2]+0.2
        du[2] =u[1]-u[2]+1.2 =#
        #= du[1] = 0.01*u[2]
        du[2] =-100.0*u[1]-100.0*u[2]+2020.0 =#

        du[1] = -4.0*u[1]+0.5*u[2]+1.0
        du[2] =u[1]-0.25*u[2]+0.25
    
    end
    tspan = (0.0,50.0)
    u0 = [100.0,0.0]
    prob = ODEProblem(funcName,u0,tspan)
    sol = solve(prob,BS3(),saveat=0.5,abstol = 1e-9, reltol = 1e-5)
    #sol = solve(prob,Rosenbrock23(),abstol = 1e-6, reltol = 1e-3)
  
    u1=(-sqrt(257)-15)/8
    u2=(sqrt(257)-15)/8
    λ1=(-sqrt(257)-17)/8
    λ2=(sqrt(257)-17)/8
    c2=(397-2*(sqrt(257)+15))/sqrt(257)
    c1=-4-c2
    b1=1.0
    b2=0.25
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+0.5*b1+b2
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+2*b1+8*b2
   #= sumTrueSqr1=0.0
   sumDiffSqr1=0.0
   relerror1=0.0
   numPoints=length(sol.t)
   for i = 1:numPoints
      ft=x2(sol.t[i])
      sumDiffSqr1+=(sol.u[i][2]-ft)*(sol.u[i][2]-ft)
      sumTrueSqr1+=ft*ft
    end
    relerror1=sqrt(sumDiffSqr1/sumTrueSqr1)
    display(relerror1) =#

    numPoints=length(sol.t)
   
    
      temp = []
      for i = 1:numPoints #each point is a taylor
        ft=x2(sol.t[i])
        push!(temp, abs(sol.u[i][2]-ft)/ft)
      end
     display(plot!(sol.t, temp,label="BS3")) 
     println("done")
     readline()

   # display(plot!(sol))
 #=    display(plot!(sol,xlims=(100,160),ylims=(-0.000002,0.000002)))
  #  display(plot!(sol,xlims=(10.0,10.0005),ylims=(-0.49971,-0.49960)))
    println("done")
    readline()  =#
end
#@btime 

odeDiffEquPackage()   =#









#= 
u1=(-sqrt(257)-15)/8
u2=(sqrt(257)-15)/8
λ1=(-sqrt(257)-17)/8
λ2=(sqrt(257)-17)/8
c2=(397-2*(sqrt(257)+15))/sqrt(257)
c1=-4-c2
b1=1.0
b2=0.25
x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+0.5*b1+b2
x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+2*b1+8*b2
#= display(plot!(x1,title="against liqss3: quan-5  ",label="true x1",xlims=(0,1),ylims=(-0.000002,12.000002)))
display(plot!(x2,label="true x2",xlims=(0,1),ylims=(-0.000002,12.000002))) =#
display(plot!(x1,title="against liqss3: quan-5  ",label="true x1",xlims=(0,20)))
display(plot!(x2,label="true x2",xlims=(0,20)))
println("press enter to exit")
readline() =#