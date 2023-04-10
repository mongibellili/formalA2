using formalqssA
#using BenchmarkTools
#using Plots
#using OrdinaryDiffEq

#include("D://QS_Solver//qss//src//models//classicProblem.jl") 
#D:\qssv01\qssv0.1\qssv01\src\models\classicProblem.jl

function test()
    odeprob = @NLodeProblemLoop begin
        u = [1.0, -1.3,0.5]
        discrete = [0.0]
        du[1] = -u[1]+u[2]
        du[3]=u[2]-u[3]
        du[2] =u[1]-u[2]   
       #=  for i in 2:3  
            du[i]=u[i-1]-u[i];
       end  =#
       #=  if 0.5-u[1]>0   #5*discrte gave error
            discrete[1]=1.0   #discrete=0.0-->type Symbol has no field args...find to personalize error msg            
        else
            discrete[1]=0.0                                    
        end =#
        
    end

    solmliqss2=QSS_Solve(odeprob,10.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
   # solliqss2=QSS_Solve_from_model(twoVarSys14,odeprob,10.0,qss3(),saveat(0.1),0.0,1e-6,1e-3)

 save_Sol(solmliqss2)
# save_Sol(solliqss3)
end
#@btime 
test()



