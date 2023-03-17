using formalqssAB
using BSON
include("/home/unknown/formalqssAB/formalqssAB/src/models/classicProblem.jl")# where you saved the model
function test1()
  odeprob = @NLodeProblem begin
        
        #----------twoVarSys1-------------
        
        u = [1.0, 0.0]
        discrete = [0.0]
        du[1] = u[2]
        du[2] =-u[1]-u[2] 
       
          
  end
  solmliqss2=QSS_Solve_from_model(twoVarSys1,odeprob,1000.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
  BSON.@save "formalAB_sys1_mliqss2_errorE-3_step0.1_ft1000.bson" solmliqss2
  solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
  BSON.@save "formalAB_sys1_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
end 
test1()
function test2()
  odeprob = @NLodeProblem begin     
        #------------twoVarSys12-----------   
        u = [1.0, 0.0]
        discrete = [0.0]
        du[1] = 0.01*u[2]
        du[2] =-100.0*u[1]-100.0*u[2]+2020.0 
  end
  solmliqss2=QSS_Solve_from_model(twoVarSys12,odeprob,1000.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
  BSON.@save "formalAB_sys2_mliqss2_errorE-3_step0.1_ft1000.bson" solmliqss2
  solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
  BSON.@save "formalAB_sys2_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
end 
test2()

#= function test3()
  odeprob = @NLodeProblem begin
        #------------twoVarSys13----------
        u = [-4.0, 4.0]
        discrete = [0.0]
        du[1] = -u[1]-u[2]+0.2
        du[2] =u[1]-u[2]+1.2    
  end
  solmliqss2=QSS_Solve_from_model(twoVarSys13,odeprob,1000.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
  BSON.@save "formalAB_sys13_mliqss2_errorE-3_step0.1_ft1000.bson" solmliqss2
  solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
  BSON.@save "formalAB_sys13_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
end 
test3() =#

function test4()
  odeprob = @NLodeProblem begin 
        #------------twoVarSys14----------
        u = [-4.0, 4.0]
        discrete = [0.0]
        du[1] = -u[1]-10.0*u[2]+0.2
        du[2] =10.0*u[1]-u[2]+1.2   
  end
  solmliqss2=QSS_Solve_from_model(twoVarSys14,odeprob,1000.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
  BSON.@save "formalAB_sys4_mliqss2_errorE-3_step0.1_ft1000.bson" solmliqss2
  solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
  BSON.@save "formalAB_sys4_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
end 
test4()

function test5()
  odeprob = @NLodeProblem begin
        #------------twoVarSys15----------
        u = [100.0, 0.0]
        discrete = [0.0]
        du[1] = -4.0*u[1]+0.5*u[2]+1.0
        du[2] =u[1]-0.25*u[2]+0.25   
  end
  solmliqss2=QSS_Solve_from_model(twoVarSys15,odeprob,1000.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
  BSON.@save "formalAB_sys5_mliqss2_errorE-3_step0.1_ft1000.bson" solmliqss2
  solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
  BSON.@save "formalAB_sys5_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
end 
test5()
function test6()
  odeprob = @NLodeProblem begin      
        #------------twoVarSys16----------
        u = [1.0, 2.3]
        discrete = [0.0]
        du[1] = -1.01*u[1]-100.0*u[2]+1000.2
        du[2] =10.1*u[1]-1.001*u[2]+0.2    
  end
  solmliqss2=QSS_Solve_from_model(twoVarSys16,odeprob,1000.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
  BSON.@save "formalAB_sys6_mliqss2_errorE-3_step0.1_ft1000.bson" solmliqss2
  solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
  BSON.@save "formalAB_sys6_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
end 
test6()