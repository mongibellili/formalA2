using relaxedqssAC


#= 

function test11()
    odeprob = @NLodeProblem begin
         #sys b11
        u = [-1.0, -2.0]
        discrete = [0.0]
        du[1] = 0.15*u[1]+0.25*u[2]
        du[2] =-7.0*u[1]-4.0*u[2]
       
    end
  
    save_prob_to_model(odeprob,"/home/unknown/models/BProblem.jl","B12") #any location you want
    #= solmliqss1=QSS_Solve(odeprob,800.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
    plotSol(solmliqss1) =#
   
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
  
    save_prob_to_model(odeprob,"/home/unknown/models/BProblem.jl","B12") #any location you want
    #= solmliqss1=QSS_Solve(odeprob,800.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
    plotSol(solmliqss1) =#
   
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
  
    save_prob_to_model(odeprob,"/home/unknown/models/BProblem.jl","B13") #any location you want
    #= solmliqss1=QSS_Solve(odeprob,800.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
    plotSol(solmliqss1) =#
   
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
  
    save_prob_to_model(odeprob,"/home/unknown/models/BProblem.jl","B14") #any location you want
   #=  solmliqss1=QSS_Solve(odeprob,800.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
    plotSol(solmliqss1) =#
   
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
  
    save_prob_to_model(odeprob,"/home/unknown/models/BProblem.jl","B15") #any location you want
    #= solmliqss1=QSS_Solve(odeprob,800.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
    plotSol(solmliqss1) =#
   
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
  
    save_prob_to_model(odeprob,"/home/unknown/models/BProblem.jl","B21") #any location you want
    #= solmliqss1=QSS_Solve(odeprob,800.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
    plotSol(solmliqss1) =#
   
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
  
    save_prob_to_model(odeprob,"/home/unknown/models/BProblem.jl","B222") #any location you want
    #= solmliqss1=QSS_Solve(odeprob,800.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
    plotSol(solmliqss1) =#
   
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
  
    save_prob_to_model(odeprob,"/home/unknown/models/BProblem.jl","B23") #any location you want
   #=  solmliqss1=QSS_Solve(odeprob,800.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
    plotSol(solmliqss1) =#
   
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
  
    save_prob_to_model(odeprob,"/home/unknown/models/BProblem.jl","B24") #any location you want
    #= solmliqss1=QSS_Solve(odeprob,800.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
    plotSol(solmliqss1)
    =#
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
  
    save_prob_to_model(odeprob,"/home/unknown/models/BProblem.jl","B25") #any location you want
   #=  solmliqss1=QSS_Solve(odeprob,800.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
    plotSol(solmliqss1) =#
   
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
  
    save_prob_to_model(odeprob,"/home/unknown/models/BProblem.jl","B31") #any location you want
   #=  solmliqss1=QSS_Solve(odeprob,800.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
    plotSol(solmliqss1) =#
   
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
  
    save_prob_to_model(odeprob,"/home/unknown/models/BProblem.jl","B32") #any location you want
   #=  solmliqss1=QSS_Solve(odeprob,800.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
    plotSol(solmliqss1) =#
   
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
  
    save_prob_to_model(odeprob,"/home/unknown/models/BProblem.jl","B33") #any location you want
   #=  solmliqss1=QSS_Solve(odeprob,800.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
    plotSol(solmliqss1) =#
   
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
  
    save_prob_to_model(odeprob,"/home/unknown/models/BProblem.jl","B34") #any location you want
   #=  solmliqss1=QSS_Solve(odeprob,800.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
    plotSol(solmliqss1) =#
   
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
    save_prob_to_model(odeprob,"/home/unknown/models/BProblem.jl","B35") #any location you want
   #=  solmliqss1=QSS_Solve(odeprob,800.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
    plotSol(solmliqss1) =#
   
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
    save_prob_to_model(odeprob,"/home/unknown/models/BProblem.jl","B41") #any location you want   
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
    save_prob_to_model(odeprob,"/home/unknown/models/BProblem.jl","B42") #any location you want
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
    save_prob_to_model(odeprob,"/home/unknown/models/BProblem.jl","B43") #any location you want 
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
  
    save_prob_to_model(odeprob,"/home/unknown/models/BProblem.jl","B44") #any location you want  
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
    save_prob_to_model(odeprob,"/home/unknown/models/BProblem.jl","B455") #any location you want 
end
test45() 
#= 
function test51()
    odeprob = @NLodeProblem begin
         #sys b51
         u = [-1.0, -2.0]
         discrete = [0.0]
         du[1] = -20.0*u[1]-80.0*u[2]
         du[2] =1.24*u[1]-0.01*u[2]
       
    end  
    save_prob_to_model(odeprob,"/home/unknown/models/BProblem.jl","B51") #any location you want   
end
test51()
 function test52()
    odeprob = @NLodeProblem begin
         #sys b52
        u = [-1.0, -2.0]
        discrete = [0.0]
        du[1] = -20.0*u[1]-80.0*u[2]+80.0
        du[2] =1.24*u[1]-0.01*u[2]+0.01      
    end
    save_prob_to_model(odeprob,"/home/unknown/models/BProblem.jl","B52") #any location you want
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
    save_prob_to_model(odeprob,"/home/unknown/models/BProblem.jl","B53") #any location you want 
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
  
    save_prob_to_model(odeprob,"/home/unknown/models/BProblem.jl","B54") #any location you want  
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
    save_prob_to_model(odeprob,"/home/unknown/models/BProblem.jl","B55") #any location you want 
end
test55()  =#