macro test(ex)
    Base.remove_linenums!(ex)
    dump(ex)
end
#= function test1()
@test begin
    u[1,2]
    discrete=[0.0]
    du[1]=u[2]-u[1]
    du[2]=u[3]-u[2]
    if u[1]>0
        discrete[1]=1.0
    end
end
end =#



function test2()
    odeprob = @test begin
        parameter1=3000.0# cache can be dynamic....parameters take this feature
        parameter2=0.00001
        u = [10.0,0.0]
        discrete = [0.0]
        du[1] =u[2]
        du[2] =-9.8#-(discrete[1])*(parameter1*u[2])
        if u[1]>0   #5*discrte gave error
            discrete[1]=0.0   #discrete=0.0-->type Symbol has no field args...find to personalize error msg            
        else
            discrete[1]=1.0                                    
        end
    end
   #=  sol= QSS_Solve(odeprob,1.0,qss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(sol)
    
    sol= QSS_Solve(odeprob,1.0,qss3(),saveat(0.01),0.0,1e-6,1e-3)
      save_Sol(sol) =#
end
test2()