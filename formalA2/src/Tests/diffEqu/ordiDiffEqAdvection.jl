using OrdinaryDiffEq
#using DifferentialEquations
#using ODEInterfaceDiffEq
#using BSON
#using qssv01
#using BenchmarkTools
using Plots
#using TimerOutputs
#reset_timer!()
function odeDiffEquPackage()
 
    function funcName(du,u,p,t)# api requires four args
   
      #= du[1] = -120.0*u[1]+10.0*u[2]+1000.0*u[1]*u[1]*(1.0-u[1]) 
       
      for i in 2:9  
          du[i]=10.0*u[i+1]+110.0*u[i-1]-120.0*u[i]+1000.0*u[i]*u[i]*(1.0-u[i]) ;
     end 

      du[10]=120.0*u[9]-120.0*u[10]+1000.0*u[10]*u[10]*(1.0-u[10])  =#
      
      #= du[1] = -120.0*u[1]+10.0*u[2]+1000.0*u[1]*u[1]*(1.0-u[1]) 
         
      for i in 2:3  
          du[i]=10.0*u[i+1]+120.0*u[i-1]-120.0*u[i]+1000.0*u[i]*u[i]*(1.0-u[i]) ;
     end 

      du[4]=120.0*u[3]-120.0*u[4]+1000.0*u[4]*u[4]*(1.0-u[4]) 
      
  
    end
    tspan = (0.0,10.0)
    #u0 =[1.0, 1.0, 1.0, 0.0, 0.0,0.0,0.0,0.0,0.0,0.0]
    u0 =[1.0, 0.0, 0.0,0.0]
    prob = ODEProblem(funcName,u0,tspan) =#

    
        du[1] = -120.0*u[1]+10.0*u[2]+1000.0*u[1]*u[1]*(1.0-u[1]) 
         
        for i in 2:99  
            du[i]=10.0*u[i+1]+110.0*u[i-1]-120.0*u[i]+1000.0*u[i]*u[i]*(1.0-u[i]) ;
       end 

        du[100]=120.0*u[99]-120.0*u[100]+1000.0*u[100]*u[100]*(1.0-u[100]) 
        
    

    
    end
    tspan = (0.0,10.0)
    u0= [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    prob = ODEProblem(funcName,u0,tspan)


    #@timeit "bs3" 
   # sol = solve(prob,BS3(),saveat=0.05,abstol = 1e-6, reltol = 1e-3)
   # @timeit "rodas5" 
    #solRodas5P = solve(prob,Rodas5P(),saveat=0.01,abstol = 1e-9, reltol = 1e-5)
    #solradau5 = solve(prob,radau5(),saveat=0.01,abstol = 1e-12, reltol = 1e-7)
     sol = solve(prob,Feagin14(),saveat=0.01,abstol = 1e-12, reltol = 1e-7)
    #=
    #display(solRodas5P);println()
     # display(solRodas5P);println()
     solFeagin14VectorN100=solFeagin14.u
    # display(roads5VectorN10)
    # display(solRodas5P(0.66)[1])
    BSON.@save "solVectAdvection_N100_Feagin14e-12.bson" solFeagin14VectorN100 =#
    #sol = solve(prob,Rosenbrock23(),abstol = 1e-6, reltol = 1e-3)
  
   
   # sol = solve(prob,BS3(),saveat=0.01,abstol = 1e-5, reltol = 1e-5)
    
    #display(solRodas5P);println()
     # display(solRodas5P);println()
    # solbs3VectorN100=solbs3.u
    # display(roads5VectorN10)
    # display(solRodas5P(0.66)[1])
    #BSON.@save "solVectAdvection_N100_bs3e-12.bson" solbs3VectorN100
   


   p1=plot!(sol,marker=(:circle),markersize=2,xlims=(0.0,10.0),ylims=(0.99995,1.00002),legend = false)
   # display(plot!(sol,xlims=(100,160),ylims=(-0.000002,0.000002)))
  #  display(plot!(sol,xlims=(10.0,10.0005),ylims=(-0.49971,-0.49960)))
  savefig(p1, "plot_adrN100_bs3_1e-5.png")
 #print_timer()
end
#@time 
odeDiffEquPackage()  



#= function funcName(du,u,p,t)# api requires four args
   
      du[1] = -120.0*u[1]+10.0*u[2]+1000.0*u[1]*u[1]*(1.0-u[1]) 
       
      for i in 2:9  
          du[i]=10.0*u[i+1]+110.0*u[i-1]-120.0*u[i]+1000.0*u[i]*u[i]*(1.0-u[i]) ;
     end 

      du[10]=120.0*u[9]-120.0*u[10]+1000.0*u[10]*u[10]*(1.0-u[10]) 
      
    
    end
    tspan = (0.0,10.0)
    u0 =[1.0, 1.0, 1.0, 0.0, 0.0,0.0,0.0,0.0,0.0,0.0]
    prob = ODEProblem(funcName,u0,tspan)
    #@timeit "bs3" 
   # sol = solve(prob,BS3(),saveat=0.05,abstol = 1e-6, reltol = 1e-3)
   # @timeit "rodas5" 
    solRodas5P = solve(prob,Rodas5P(),saveat=0.01,abstol = 1e-9, reltol = 1e-5)
    display(solRodas5P)
    BSON.@save "solRodas5Pe-9.bson" solRodas5P =#

