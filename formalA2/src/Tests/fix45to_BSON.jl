using formalqss
using BSON
include("/home/unknown/formalA2/formalA2/src/models/InvestigateBProblem.jl") 
 #=function test11()
    odeprob = @NLodeProblem begin
         #sys b11
        u = [-1.0, -2.0]
        discrete = [0.0]
        du[1] = 0.15*u[1]+0.25*u[2]
        du[2] =-7.0*u[1]-4.0*u[2]
       
    end
   #=  solmliqss2=QSS_Solve_from_model(B11,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(solmliqss2,"intgr: noDir mquantizer:hyb_ft1000") =#
   # save_SolZoomed(solmliqss2,"zoomed_intgr: noDir mquantizer:dir_ft1000")
   # save_prob_to_model(odeprob,"D://QS_Solver//qss//src//models//InvestigateBProblem.jl","B12") #any location you want
    solmliqss2=QSS_Solve_from_model(B11,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(solmliqss2)
    BSON.@save "sysB/formal_sysB11_mliqss2_errorE-3_step0.01_ft1000.bson" solmliqss2
    solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
    BSON.@save "sysB/formal_sysB11_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
   
end
test11() 
 function test12()
    odeprob = @NLodeProblem begin
         #sys b11
        u = [-1.0, -2.0]
        discrete = [0.0]
        du[1] = 0.15*u[1]+0.25*u[2]-0.25
        du[2] =-7.0*u[1]-4.0*u[2]+4.0
       
    end
  
   #=  solmliqss2=QSS_Solve_from_model(B51,odeprob,12.5,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
 #  save_Sol(solmliqss2,"intgr: noDir mquantizer:dir_ft1.5")
   save_SolZoomed(solmliqss2,"zoomed_intgr: noDir mquantizer:dir_ft12.5") =#
   #= solmliqss2=QSS_Solve_from_model(B12,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
   save_Sol(solmliqss2,"intgr: noDir mquantizer:hyb_ft1000") =#
   #save_SolZoomed(solmliqss2,"zoomed_intgr: noDir mquantizer:dir_ft100")


   # save_prob_to_model(odeprob,"D://QS_Solver//qss//src//models//InvestigateBProblem.jl","B12") #any location you want
   solmliqss2=QSS_Solve_from_model(B12,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
   save_Sol(solmliqss2)
   BSON.@save "sysB/formal_sysB12_mliqss2_errorE-3_step0.01_ft1000.bson" solmliqss2
   solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
   BSON.@save "sysB/formal_sysB12_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
   
end
test12() 

 function test13()
    odeprob = @NLodeProblem begin
         #sys b13
        u = [-1.0, -2.0]
        discrete = [0.0]
        du[1] = 0.15*u[1]+0.25*u[2]-5.0
        du[2] =-7.0*u[1]-4.0*u[2]+80.0
       
       

    end
    #=  solmliqss2=QSS_Solve_from_model(B51,odeprob,12.5,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
 #  save_Sol(solmliqss2,"intgr: noDir mquantizer:dir_ft1.5")
   save_SolZoomed(solmliqss2,"zoomed_intgr: noDir mquantizer:dir_ft12.5") =#
  #=  solmliqss2=QSS_Solve_from_model(B13,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
   save_Sol(solmliqss2,"intgr: noDir mquantizer:Hyb_ft1000") =#
  # save_SolZoomed(solmliqss2,"zoomed_intgr: noDir mquantizer:dir_ft30")



    solmliqss2=QSS_Solve_from_model(B13,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
   save_Sol(solmliqss2)
   BSON.@save "sysB/formal_sysB13_mliqss2_errorE-3_step0.01_ft1000.bson" solmliqss2
   solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
   BSON.@save "sysB/formal_sysB13_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
   
end
test13()  

 function test14()
    odeprob = @NLodeProblem begin
         #sys b14
        u = [-1.0, -2.0]
        discrete = [0.0]
        du[1] = 0.15*u[1]+0.25*u[2]-0.0025
        du[2] =-7.0*u[1]-4.0*u[2]+0.04
        
       

    end
    #= solmliqss2=QSS_Solve_from_model(B14,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
   # save_Sol(solmliqss2,"intgr: noDir mquantizer:dir_ft1000")
    save_SolZoomed(solmliqss2,"zoomed_intgr: noDir mquantizer:hyb_ft1000") =#
    solmliqss2=QSS_Solve_from_model(B14,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
   save_Sol(solmliqss2)
   BSON.@save "sysB/formal_sysB14_mliqss2_errorE-3_step0.01_ft1000.bson" solmliqss2
   solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
   BSON.@save "sysB/formal_sysB14_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
   
end

test14()

function test15()
    odeprob = @NLodeProblem begin
         
        #sys b15
        u = [-1.0, -2.0]
        discrete = [0.0]
        du[1] = 0.15*u[1]+0.25*u[2]+0.05
        du[2] =-7.0*u[1]-4.0*u[2]-10.0
       

    end
    #= solmliqss2=QSS_Solve_from_model(B15,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(solmliqss2,"intgr: noDir mquantizer:hyb_ft1000") =#
  #  save_SolZoomed(solmliqss2,"zoomed_intgr: noDir mquantizer:dir_ft1000")
    solmliqss2=QSS_Solve_from_model(B15,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
   save_Sol(solmliqss2)
   BSON.@save "sysB/formal_sysB15_mliqss2_errorE-3_step0.01_ft1000.bson" solmliqss2
   solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
   BSON.@save "sysB/formal_sysB15_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
   
end
test15()
function test21()
    odeprob = @NLodeProblem begin
        
         #sys b21
        u = [-1.0, -2.0]
        discrete = [0.0]
        du[1] = 0.1*u[1]+0.25*u[2]
        du[2] =-4.0*u[1]-7.0*u[2]


    end
  #=   solmliqss2=QSS_Solve_from_model(B21,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(solmliqss2,"intgr: noDir mquantizer:hyb_ft1000") =#
  #  save_SolZoomed(solmliqss2,"zoomed_intgr: noDir mquantizer:dir_ft1000")
    solmliqss2=QSS_Solve_from_model(B21,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(solmliqss2)
    BSON.@save "sysB/formal_sysB21_mliqss2_errorE-3_step0.01_ft1000.bson" solmliqss2
    solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
    BSON.@save "sysB/formal_sysB21_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp

   
end
test21()
 function test22()
    odeprob = @NLodeProblem begin
        
         #sys b22
        u = [-1.0, -2.0]
        discrete = [0.0]
        du[1] = 0.1*u[1]+0.25*u[2]-0.25
        du[2] =-4.0*u[1]-7.0*u[2]+7.0

       

    end
   #=  solmliqss2=QSS_Solve_from_model(B22,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(solmliqss2,"intgr: noDir mquantizer:hyb_ft1000") =#
  #  save_SolZoomed(solmliqss2,"zoomed_intgr: noDir mquantizer:dir_ft1000")
    solmliqss2=QSS_Solve_from_model(B22,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(solmliqss2)
    BSON.@save "sysB/formal_sysB22_mliqss2_errorE-3_step0.01_ft1000.bson" solmliqss2
    solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
    BSON.@save "sysB/formal_sysB22_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
   
end
test22() 

function test23()
    odeprob = @NLodeProblem begin
       
         #sys b25
        u = [-1.0, -2.0]
        discrete = [0.0]
        du[1] = 0.1*u[1]+0.25*u[2]-5.0
        du[2] =-4.0*u[1]-7.0*u[2]+140.0

       
       

    end
    #= solmliqss2=QSS_Solve_from_model(B23,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(solmliqss2,"intgr: noDir mquantizer:hyb_ft1000") =#
   # save_SolZoomed(solmliqss2,"zoomed_intgr: noDir mquantizer:dir_ft1000")
    solmliqss2=QSS_Solve_from_model(B23,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(solmliqss2)
    BSON.@save "sysB/formal_sysB23_mliqss2_errorE-3_step0.01_ft1000.bson" solmliqss2
    solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
    BSON.@save "sysB/formal_sysB23_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
   
end
test23()
function test24()
    odeprob = @NLodeProblem begin
        
         #sys b24
        u = [-1.0, -2.0]
        discrete = [0.0]
        du[1] = 0.1*u[1]+0.25*u[2]-0.0025
        du[2] =-4.0*u[1]-7.0*u[2]+0.07

       

    end
   #=  solmliqss2=QSS_Solve_from_model(B24,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(solmliqss2,"intgr: noDir mquantizer:hyb_ft1000") =#
   # save_SolZoomed(solmliqss2,"zoomed_intgr: noDir mquantizer:dir_ft1000")
    solmliqss2=QSS_Solve_from_model(B24,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(solmliqss2)
    BSON.@save "sysB/formal_sysB24_mliqss2_errorE-3_step0.01_ft1000.bson" solmliqss2
    solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
    BSON.@save "sysB/formal_sysB24_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
end
test24()
function test25()
    odeprob = @NLodeProblem begin
         
         #sys b25
        u = [-1.0, -2.0]
        discrete = [0.0]
        du[1] = 0.1*u[1]+0.25*u[2]-0.05
        du[2] =-4.0*u[1]-7.0*u[2]-1.0

      
       

    end
   #=  solmliqss2=QSS_Solve_from_model(B25,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(solmliqss2,"intgr: noDir mquantizer:hyb_ft1000") =#
   # save_SolZoomed(solmliqss2,"zoomed_intgr: noDir mquantizer:dir_ft1000")
    solmliqss2=QSS_Solve_from_model(B25,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
   save_Sol(solmliqss2)
   BSON.@save "sysB/formal_sysB25_mliqss2_errorE-3_step0.01_ft1000.bson" solmliqss2
   solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
   BSON.@save "sysB/formal_sysB25_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
   
end
test25()
function test31()
    odeprob = @NLodeProblem begin
         
         #sys b31
        u = [-1.0, -2.0]
        discrete = [0.0]
        du[1] = -0.25*u[1]+0.15*u[2]
        du[2] =-4.0*u[1]-7.0*u[2]

      
       

    end

  #=  solmliqss2=QSS_Solve_from_model(B31,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(solmliqss2,"intgr: noDir mquantizer:hyb_ft1000") =#
  #  save_SolZoomed(solmliqss2,"zoomed_intgr: noDir mquantizer:dir_ft1000")
    solmliqss2=QSS_Solve_from_model(B31,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
   save_Sol(solmliqss2)
   BSON.@save "sysB/formal_sysB31_mliqss2_errorE-3_step0.01_ft1000.bson" solmliqss2
   solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
   BSON.@save "sysB/formal_sysB31_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
   
end
test31()

function test32()
    odeprob = @NLodeProblem begin
         
         #sys b31
        u = [-1.0, -2.0]
        discrete = [0.0]
        du[1] = -0.25*u[1]+0.15*u[2]-0.15
        du[2] =-4.0*u[1]-7.0*u[2]+7.0

      
       

    end
    #= solmliqss2=QSS_Solve_from_model(B32,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(solmliqss2,"intgr: noDir mquantizer:hyb_ft1000") =#
   # save_SolZoomed(solmliqss2,"zoomed_intgr: noDir mquantizer:ddx_ft100")
    solmliqss2=QSS_Solve_from_model(B32,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(solmliqss2)
    BSON.@save "sysB/formal_sysB32_mliqss2_errorE-3_step0.01_ft1000.bson" solmliqss2
    solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
    BSON.@save "sysB/formal_sysB32_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
   
end
test32()
function test33()
    odeprob = @NLodeProblem begin
         
         #sys b31
        u = [-1.0, -2.0]
        discrete = [0.0]
        du[1] = -0.25*u[1]+0.15*u[2]-3.0
        du[2] =-4.0*u[1]-7.0*u[2]+140.0

      
       

    end
  #=   solmliqss2=QSS_Solve_from_model(B33,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(solmliqss2,"intgr: noDir mquantizer:hyb_ft1000") =#
    #save_SolZoomed(solmliqss2,"zoomed_intgr: noDir mquantizer:ddx_ft100")
    solmliqss2=QSS_Solve_from_model(B33,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
   save_Sol(solmliqss2)
   BSON.@save "sysB/formal_sysB33_mliqss2_errorE-3_step0.01_ft1000.bson" solmliqss2
   solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
   BSON.@save "sysB/formal_sysB33_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
   
end
test33()

function test34()
    odeprob = @NLodeProblem begin
         
         #sys b31
        u = [-1.0, -2.0]
        discrete = [0.0]
        du[1] = -0.25*u[1]+0.15*u[2]-0.0015
        du[2] =-4.0*u[1]-7.0*u[2]+0.07

      
       

    end
    #= solmliqss2=QSS_Solve_from_model(B34,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(solmliqss2,"intgr: noDir mquantizer:hyb_ft1000") =#
    #save_SolZoomed(solmliqss2,"zoomed_intgr: noDir mquantizer:dir_ft1000")
    solmliqss2=QSS_Solve_from_model(B34,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
   save_Sol(solmliqss2)
   BSON.@save "sysB/formal_sysB34_mliqss2_errorE-3_step0.01_ft1000.bson" solmliqss2
   solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
   BSON.@save "sysB/formal_sysB34_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
   
end
test34()
function test35()
    odeprob = @NLodeProblem begin
         
         #sys b35
        u = [-1.0, -2.0]
        discrete = [0.0]
        du[1] = -0.25*u[1]+0.15*u[2]-0.65
        du[2] =-4.0*u[1]-7.0*u[2]-1.0

      
       

    end
    #= solmliqss2=QSS_Solve_from_model(B35,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(solmliqss2,"intgr: noDir mquantizer:hyb_ft1000") =#
   # save_SolZoomed(solmliqss2,"zoomed_intgr: noDir mquantizer:dir_ft1000")
    solmliqss2=QSS_Solve_from_model(B35,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
   save_Sol(solmliqss2)
   BSON.@save "sysB/formal_sysB35_mliqss2_errorE-3_step0.01_ft1000.bson" solmliqss2
   solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
   BSON.@save "sysB/formal_sysB35_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
   
end
test35()
function test41()
    odeprob = @NLodeProblem begin
         #sys b41
        u = [-1.0, -2.0]
        discrete = [0.0]
        du[1] = -0.01*u[1]+1.24*u[2]
        du[2] =-80.0*u[1]-20.0*u[2]
       
    end  
   #=  solmliqss2=QSS_Solve_from_model(B41,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(solmliqss2,"intgr: noDir mquantizer:hyb_ft1000") =#
   # save_SolZoomed(solmliqss2,"zoomed_intgr: noDir mquantizer:dir_ft1000")
    solmliqss2=QSS_Solve_from_model(B41,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
   save_Sol(solmliqss2)
   BSON.@save "sysB/formal_sysB41_mliqss2_errorE-3_step0.01_ft1000.bson" solmliqss2
   solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
   BSON.@save "sysB/formal_sysB41_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
end
test41()

 function test42()
    odeprob = @NLodeProblem begin
         #sys b42
        u = [-1.0, -2.0]
        discrete = [0.0]
        du[1] = -0.01*u[1]+1.24*u[2]-1.24
        du[2] =-80.0*u[1]-20.0*u[2]+20.0      
    end
    #= solmliqss2=QSS_Solve_from_model(B42,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(solmliqss2,"intgr: noDir mquantizer:hyb_ft1000") =#
   # save_SolZoomed(solmliqss2,"zoomed_intgr: noDir mquantizer:ddx_ft100")
    solmliqss2=QSS_Solve_from_model(B42,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
   save_Sol(solmliqss2)
   BSON.@save "sysB/formal_sysB42_mliqss2_errorE-3_step0.01_ft1000.bson" solmliqss2
   solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
   BSON.@save "sysB/formal_sysB42_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
end
test42()
function test43()
    odeprob = @NLodeProblem begin
         #sys b43
        u = [-1.0, -2.0]
        discrete = [0.0]
        du[1] = -0.01*u[1]+1.24*u[2]-24.8
        du[2] =-80.0*u[1]-20.0*u[2]+400.0
    end 
    #= solmliqss2=QSS_Solve_from_model(B43,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(solmliqss2,"intgr: noDir mquantizer:hyb_ft1000") =#
  #  save_SolZoomed(solmliqss2,"zoomed_intgr: noDir mquantizer:ddx_ft100")
    solmliqss2=QSS_Solve_from_model(B43,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
   save_Sol(solmliqss2)
   BSON.@save "sysB/formal_sysB43_mliqss2_errorE-3_step0.01_ft1000.bson" solmliqss2
   solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
   BSON.@save "sysB/formal_sysB43_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
end
test43()  

function test44()
    odeprob = @NLodeProblem begin
         #sys b44
        u = [-1.0, -2.0]
        discrete = [0.0]
        du[1] = -0.01*u[1]+1.24*u[2]-0.0124
        du[2] =-80.0*u[1]-20.0*u[2]+0.2
    end
    #= solmliqss2=QSS_Solve_from_model(B44,odeprob,20.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(solmliqss2,"hyb_ft10")
    save_SolZoomed(solmliqss2,"zoomed_hyb_ft100") =#
    solmliqss2=QSS_Solve_from_model(B44,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
   save_Sol(solmliqss2)
   BSON.@save "sysB/formal_sysB44_mliqss2_errorE-3_step0.01_ft1000.bson" solmliqss2
   solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
   BSON.@save "sysB/formal_sysB44_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
 
end
test44() =#
  function test45()
    odeprob = @NLodeProblem begin         
        #sys b45
        u = [-1.0, -2.0]
        discrete = [0.0]
        du[1] = -0.01*u[1]+1.24*u[2]-1.26
        du[2] =-80.0*u[1]-20.0*u[2]-140.0
    end
   #=  solmliqss2=QSS_Solve_from_model(B45,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(solmliqss2,"intgr: noDir mquantizer:hyb_ft1000") =#
  #  save_SolZoomed(solmliqss2,"zoomed_intgr: noDir mquantizer:dir_ft1000")
    solmliqss2=QSS_Solve_from_model(B45,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
   save_Sol(solmliqss2)
   BSON.@save "sysB/formal_sysB45_mliqss2_errorE-3_step0.01_ft1000.bson" solmliqss2
   solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
   BSON.@save "sysB/formal_sysB45_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
end
test45() 
#=
##########################################

 function test51()
    odeprob = @NLodeProblem begin
         #sys b51
        u = [-1.0, -2.0]
        discrete = [0.0]
        du[1] = -20.0*u[1]-80.0*u[2]
        du[2] =1.24*u[1]-0.01*u[2]
       
    end  
     #= solmliqss2=QSS_Solve_from_model(B51,odeprob,12.5,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
 #  save_Sol(solmliqss2,"intgr: noDir mquantizer:dir_ft1.5")
   save_SolZoomed(solmliqss2,"zoomed_intgr: noDir mquantizer:dir_ft12.5")  =#
   solmliqss2=QSS_Solve_from_model(B51,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
   save_Sol(solmliqss2) 
   #save_SolZoomed(solmliqss2,"zoomed_intgr: noDir mquantizer:dir_ft100")



    BSON.@save "sysB/formal_sysB51_mliqss2_errorE-3_step0.01_ft1000.bson" solmliqss2
   solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
   BSON.@save "sysB/formal_sysB51_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp

end
#@time
 test51() 
function test52()
    odeprob = @NLodeProblem begin
         #sys b52
        u = [-1.0, -2.0]
        discrete = [0.0]
        du[1] = -20.0*u[1]-80.0*u[2]+80.0
        du[2] =1.24*u[1]-0.01*u[2]+0.01     
    end
    #= solmliqss2=QSS_Solve_from_model(B52,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(solmliqss2,"intgr: noDir mquantizer:hyb_ft1000") =#
    #save_SolZoomed(solmliqss2,"zoomed_intgr: noDir mquantizer:dir_ft1000")
    solmliqss2=QSS_Solve_from_model(B52,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
   save_Sol(solmliqss2)
   BSON.@save "sysB/formal_sysB52_mliqss2_errorE-3_step0.01_ft1000.bson" solmliqss2
   solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
   BSON.@save "sysB/formal_sysB52_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
end
test52()
function test53()
    odeprob = @NLodeProblem begin
         #sys b53
        u = [-1.0, -2.0]
        discrete = [0.0]
        du[1] = -20.0*u[1]-80.0*u[2]+1600.0
        du[2] =1.24*u[1]-0.01*u[2]+0.2
    end 
    #= solmliqss2=QSS_Solve_from_model(B53,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(solmliqss2,"intgr: noDir mquantizer:hyb_ft1000") =#
   # save_SolZoomed(solmliqss2,"zoomed_intgr: noDir mquantizer:dir_ft10")
    solmliqss2=QSS_Solve_from_model(B53,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
  save_Sol(solmliqss2)
   BSON.@save "sysB/formal_sysB53_mliqss2_errorE-3_step0.01_ft1000.bson" solmliqss2
   solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
   BSON.@save "sysB/formal_sysB53_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
end
test53()

function test54()
    odeprob = @NLodeProblem begin
         #sys b54
        u = [-1.0, -2.0]
        discrete = [0.0]
        du[1] = -20.0*u[1]-80.0*u[2]+0.8
        du[2] =1.24*u[1]-0.01*u[2]+0.0001
    end
 #=    solmliqss2=QSS_Solve_from_model(B54,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(solmliqss2,"intgr: noDir mquantizer:hyb_ft1000") =#
   # save_SolZoomed(solmliqss2,"zoomed_intgr: noDir mquantizer:dir_ft1000")
    solmliqss2=QSS_Solve_from_model(B54,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
   save_Sol(solmliqss2)
   BSON.@save "sysB/formal_sysB54_mliqss2_errorE-3_step0.01_ft1000.bson" solmliqss2
   solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
   BSON.@save "sysB/formal_sysB54_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
end
test54() 
function test55()
    odeprob = @NLodeProblem begin         
        #sys b55
        u = [-1.0, -2.0]
        discrete = [0.0]
        du[1] = -20.0*u[1]-80.0*u[2]+40.0
        du[2] =1.24*u[1]-0.01*u[2]+2.49
    end
  #=   solmliqss2=QSS_Solve_from_model(B55,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(solmliqss2,"intgr: noDir mquantizer:hyb_ft1000") =#
  # save_SolZoomed(solmliqss2,"zoomed_intgr: noDir mquantizer:dir_ft1000")
    solmliqss2=QSS_Solve_from_model(B55,odeprob,1000.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
   save_Sol(solmliqss2)
   BSON.@save "sysB/formal_sysB55_mliqss2_errorE-3_step0.01_ft1000.bson" solmliqss2
   solmliqss2Interp=solInterpolated(solmliqss2,0.1,1000.0)
   BSON.@save "sysB/formal_sysB55_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
end
test55()  =#