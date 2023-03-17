using qssmodifiedv02

include("D://QS_Solver//qss//src//models//InvestigateA0Problem.jl") 
function testA011()
    odeprob = @NLodeProblem begin
        #--------------A011-----------------
        u = [-1.0, 2.0]
        discrete = [0.0]
        du[1] = 2.0*u[2]
        du[2] =-0.1*u[1]-10.0*u[2]
      
    end
  
   
    solmliqss1=QSS_Solve_from_model(A011,odeprob,1000.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
    save_Sol(solmliqss1)
   
end
testA011()
function testA012()
    odeprob = @NLodeProblem begin
       
       #--------------A012-----------------
        u = [-1.0, 2.0]
        discrete = [0.0]
        du[1] = 2.0*u[2]-2.0
        du[2] =-0.1*u[1]-10.0*u[2]+10.0 
      
    end
  
   
    solmliqss1=QSS_Solve_from_model(A012,odeprob,1000.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
    save_Sol(solmliqss1)
   
end
testA012()
function testA013()
    odeprob = @NLodeProblem begin
       
        #--------------A013-----------------
        u = [-1.0, 2.0]
        discrete = [0.0]
        du[1] = 2.0*u[2]-40.0
        du[2] =-0.1*u[1]-10.0*u[2]+200.0
       

    end
  
   
    solmliqss1=QSS_Solve_from_model(A013,odeprob,1000.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
    save_Sol(solmliqss1)
   
end
testA013()
function testA014()
    odeprob = @NLodeProblem begin
        
        #--------------A014-----------------
        u = [-1.0, 2.0]
        discrete = [0.0]
        du[1] = 2.0*u[2]-0.02
        du[2] =-0.1*u[1]-10.0*u[2]+0.1
      
       
    end
  
   
    solmliqss1=QSS_Solve_from_model(A014,odeprob,1000.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
    save_Sol(solmliqss1)
   
end
testA014()
function testA015()
    odeprob = @NLodeProblem begin
      
      
        #------------------A015--------------------
        u = [-1.0, -2.0]
        discrete = [0.0]
        du[1] = 2.0*u[2]-2.0
        du[2] =-0.1*u[1]-10.0*u[2]+9.8
       

    end
  
   
    solmliqss1=QSS_Solve_from_model(A015,odeprob,1000.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
    save_Sol(solmliqss1)
   
end
testA015()
function testA021()
    odeprob = @NLodeProblem begin
       

         #--------------A021-----------------
        u = [-1.0, 2.0]
        discrete = [0.0]
        du[1] = 0.25*u[2]
        du[2] =-8.0*u[1]-4.0*u[2]
        
    end
  
   
    solmliqss1=QSS_Solve_from_model(A021,odeprob,1000.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
    save_Sol(solmliqss1)
   
end
testA021()
function testA022()
    odeprob = @NLodeProblem begin

         #--------------A022-----------------
         u = [-1.0, 2.0]
         discrete = [0.0]
         du[1] = 0.25*u[2]-0.25
         du[2] =-8.0*u[1]-4.0*u[2]+4.0

    end
 
    solmliqss1=QSS_Solve_from_model(A022,odeprob,1000.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
    save_Sol(solmliqss1)
   
end
testA022()
function testA023()
    odeprob = @NLodeProblem begin
       
          #--------------A023-----------------
          u = [-1.0, 2.0]
          discrete = [0.0]
          du[1] = 0.25*u[2]-5.0
          du[2] =-8.0*u[1]-4.0*u[2]+80.0         
    end
   
    solmliqss1=QSS_Solve_from_model(A023,odeprob,1000.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
    save_Sol(solmliqss1)
   
end
testA023()
function testA024()
    odeprob = @NLodeProblem begin
      
          #--------------A024-----------------
          u = [-1.0, 2.0]
          discrete = [0.0]
          du[1] = 0.25*u[2]-0.0025
          du[2] =-8.0*u[1]-4.0*u[2]+0.04
     end
  
    solmliqss1=QSS_Solve_from_model(A024,odeprob,1000.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
    save_Sol(solmliqss1)
   
end
testA024()
function testA025()
    odeprob = @NLodeProblem begin
       
          #--------------A025-----------------
          u = [-1.0, 2.0]
          discrete = [0.0]
          du[1] = 0.25*u[2]-0.25
          du[2] =-8.0*u[1]-4.0*u[2]-12.0

    end
  
   
    solmliqss1=QSS_Solve_from_model(A025,odeprob,1000.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
    save_Sol(solmliqss1)
   
end
testA025()
