

using qssv01
#= using Plots;
gr(); =#
function test()
    odeprob = @NLodeProblem begin
        parameter2=30.0# why 300 bad
        parameter1=1e5
        u = [10.0,0.0]
        discrete = [0.0]
        du[1] =u[2]
        du[2] =-9.8-(discrete[1])*(parameter1*u[1]+parameter2*u[2])
        if u[1]>0   #5*discrte gave error
            discrete[1]=0.0   #discrete=0.0-->type Symbol has no field args...find to personalize error msg            
        else
            discrete[1]=1.0                                    
        end
    end
    sol= QSS_Solve(odeprob,2.0,qss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(sol)
    
   #=  sol= QSS_Solve(odeprob,1.0,qss3(),saveat(0.01),0.0,1e-6,1e-3)
      save_Sol(sol) =#
end
test()
