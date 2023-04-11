using OrdinaryDiffEq
#using DifferentialEquations
#using ODEInterfaceDiffEq
using BSON
#using qssv01
#using BenchmarkTools
using Plots
#using TimerOutputs
#reset_timer!()
function odeDiffEquPackage()
    function funcName(du,u,p,t)# api requires four args
        du[1] = -120.0*u[1]+10.0*u[2]+1000.0*u[1]*u[1]*(1.0-u[1]) 
        for i in 2:9  
            du[i]=10.0*u[i+1]+110.0*u[i-1]-120.0*u[i]+1000.0*u[i]*u[i]*(1.0-u[i]) ;
       end 
        du[10]=120.0*u[9]-120.0*u[10]+1000.0*u[10]*u[10]*(1.0-u[10]) 
    end
    tspan = (0.0,10.0)
    u0= [1.0, 1.0, 1.0, 0.0, 0.0,0.0,0.0,0.0,0.0,0.0]
    prob = ODEProblem(funcName,u0,tspan)

    #@timeit "bs3" 
   # sol = solve(prob,BS3(),saveat=0.05,abstol = 1e-6, reltol = 1e-3)
   # @timeit "rodas5" 
    #solRodas5P = solve(prob,Rodas5P(),saveat=0.01,abstol = 1e-9, reltol = 1e-5)
    #solradau5 = solve(prob,radau5(),saveat=0.01,abstol = 1e-12, reltol = 1e-7)
    solFeagin14VectorN10 = solve(prob,Feagin14(),saveat=0.01,abstol = 1e-13, reltol = 1e-7)

#= temp=[]
for i in solFeagin14VectorN10.u
   push!(temp,i[1])
end
p1=plot()
 p1=plot!(p1,solFeagin14VectorN10.t,temp,marker=(:circle),markersize=2,xlims=(0.0,10.0),ylims=(0.87585,0.87588),label="x1")
 savefig(p1, "plot_x1_adrN10_feagin14_1e-5.png") =#
#########################################################################################################
#=  temp=[]
 for i in solFeagin14VectorN10.u
    push!(temp,i[2])
 end
 p1=plot()
  p1=plot!(p1,solFeagin14VectorN10.t,temp,marker=(:circle),markersize=2,xlims=(0.0,10.0),ylims=(0.98752,0.987523),label="x2")
  savefig(p1, "plot_x2_adrN10_feagin14_1e-5.png") =#
#########################################################################################################
#=   temp=[]
  for i in solFeagin14VectorN10.u
     push!(temp,i[3])
  end
   p1=plot()
   p1=plot!(p1,solFeagin14VectorN10.t,temp,marker=(:circle),markersize=2,xlims=(0.0,10.0),ylims=(0.99875,0.99878),label="x3")
   savefig(p1, "plot_x3_adrN10_feagin14_1e-5.png") =#
   #########################################################################################################
   temp=[]
   for i in solFeagin14VectorN10.u
      push!(temp,i[4])
   end
     p1=plot()
    p1=plot!(p1,solFeagin14VectorN10.t,temp,marker=(:circle),markersize=2,xlims=(0.0,10.0),ylims=(0.99987,0.9999),label="x4")
    savefig(p1, "plot_x4_adrN10_feagin14_1e-5.png")
#########################################################################################################
  #=  temp=[]
    for i in solFeagin14VectorN10.u
       push!(temp,i[5])
    end
    
    p1=plot()
     p1=plot!(p1,solFeagin14VectorN10.t,temp,marker=(:circle),markersize=2,xlims=(0.0,10.0),ylims=(0.9999855,0.99999),label="x5")
     savefig(p1, "plot_x5_adrN10_feagin14_1e-5.png")
#########################################################################################################
     temp=[]
     for i in solFeagin14VectorN10.u
        push!(temp,i[6])
     end
     p1=plot()
      p1=plot!(p1,solFeagin14VectorN10.t,temp,marker=(:circle),markersize=2,xlims=(0.0,10.0),ylims=(0.999998,1.0),label="x6")
      savefig(p1, "plot_x6_adrN10_feagin14_1e-5.png")

 =#

 #########################################################################################################
#=  temp=[]
 for i in solFeagin14VectorN10.u
    push!(temp,i[7])
 end
 p1=plot()
  p1=plot!(p1,solFeagin14VectorN10.t,temp,marker=(:circle),markersize=2,xlims=(0.0,10.0),ylims=(0.99999,1.00001),label="x7")
  savefig(p1, "plot_x7_adrN10_feagin14_1e-5.png")

  #########################################################################################################
  temp=[]
  for i in solFeagin14VectorN10.u
     push!(temp,i[8])
  end
  p1=plot()
   p1=plot!(p1,solFeagin14VectorN10.t,temp,marker=(:circle),markersize=2,xlims=(0.0,10.0),ylims=(0.99999,1.00001),label="x8")
   savefig(p1, "plot_x8_adrN10_feagin14_1e-5.png")

   #########################################################################################################
   temp=[]
   for i in solFeagin14VectorN10.u
      push!(temp,i[9])
   end
   p1=plot()
    p1=plot!(p1,solFeagin14VectorN10.t,temp,marker=(:circle),markersize=2,xlims=(0.0,10.0),ylims=(0.99999,1.00001),label="x9")
    savefig(p1, "plot_x9_adrN10_feagin14_1e-5.png")


#########################################################################################################
temp=[]
for i in solFeagin14VectorN10.u
   push!(temp,i[10])
end
p1=plot()
 p1=plot!(p1,solFeagin14VectorN10.t,temp,marker=(:circle),markersize=2,xlims=(0.0,10.0),ylims=(0.99999,1.00001),label="x10")
 savefig(p1, "plot_x10_adrN10_feagin14_1e-5.png")
 =#







 #=   p1=plot!(solFeagin14VectorN10,marker=(:circle),markersize=2,xlims=(0.0,10.0),ylims=(0.9999995,1.0000002))
  savefig(p1, "plot_adrN10_feagin14_1e-5.png") =#

end
#@time 
odeDiffEquPackage()  


