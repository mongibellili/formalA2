
using XLSX
using formalqss
using StaticArrays
using BSON

#include("/home/unknown/formalqssAB/formalqssAB/src/models/InvestigateBProblem.jl") 
#=  function test11()
    BSON.@load "sysB/formal_sysB11_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
   
   
    u1, u2 = -0.06805253887502617, -0.5248046039821167
    λ1, λ2 = -3.523632227874817, -0.32636777212518275
   
    c1, c2 = -4.487356192869495, 2.4873561928694947
    xp1, xp2 = 0.0, 0.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
   #=  saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
   #=  er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)  =# 

     ("B11",er1,er2)
   
end
err11=test11()
function test12()
    BSON.@load "sysB/formal_sysB12_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    u1, u2 = -0.06805253887502617, -0.5248046039821167
    λ1, λ2 = -3.523632227874817, -0.32636777212518275
   
    c1, c2 = -5.636348488851935, 2.6363484888519344
    xp1, xp2 = 1.1102230246251565e-16, 0.9999999999999999
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
    #= saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
    #= er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)   =#

     ("B12",er1,er2)
   
end
err12=test12()
function test13()
    BSON.@load "sysB/formal_sysB13_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    u1, u2 = -0.06805253887502617, -0.5248046039821167
    λ1, λ2 = -3.523632227874817, -0.32636777212518275
   
    c1, c2 = -27.467202112518294, 5.467202112518295
    xp1, xp2 = 3.552713678800501e-15, 20.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
   #=  saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
    #= er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)   =#

     ("B13",er1,er2)
   
end
err13=test13()
function test14()
    BSON.@load "sysB/formal_sysB14_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    u1, u2 = -0.06805253887502617, -0.5248046039821167
    λ1, λ2 = -3.523632227874817, -0.32636777212518275
   
    c1, c2 = -4.498846115829319, 2.488846115829319
    xp1, xp2 = 0.0, 0.009999999999999998
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
   #=  saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
  #=   er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)  =# 

     ("B14",er1,er2)
end
err14=test14()
function test15()
    BSON.@load "sysB/formal_sysB15_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    u1, u2 = -0.06805253887502617, -0.5248046039821167
    λ1, λ2 = -3.523632227874817, -0.32636777212518275
   
    c1, c2 = -1.2576052870427024, -1.7423947129572968
    xp1, xp2 = -2.0, 0.9999999999999993
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
   #=  saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
   #=  er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)   =#

     ("B15",er1,er2)
   
end
err15=test15()

function test21()
    BSON.@load "sysB/formal_sysB21_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    u1, u2 = -0.03593893348744506, -1.739061066512555
    λ1, λ2 = -6.856244266050219, -0.04375573394978005
   
    c1, c2 = -2.6293605409677845, 0.6293605409677845
    xp1, xp2 = 0.0, 0.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
  #=   saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
    #= er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,60.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,60.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,60.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,60.0,1000.0)   =#

     ("B21",er1,er2)
   
end
err21=test21()
function test22()
    BSON.@load "sysB/formal_sysB22_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    u1, u2 = -0.03593893348744506, -1.739061066512555
    λ1, λ2 = -6.856244266050219, -0.04375573394978005
   
    c1, c2 = -3.6504623356016257, 0.6504623356016255
    xp1, xp2 = 0.0, 1.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
   #=  saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
    #= er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,60.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,60.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,60.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,60.0,1000.0)   =#

     ("B22",er1,er2)
   
end
err22=test22()
function test23()
    BSON.@load "sysB/formal_sysB23_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    u1, u2 = -0.03593893348744506, -1.739061066512555
    λ1, λ2 = -6.856244266050219, -0.04375573394978005
   
    c1, c2 = -23.0513964336446, 1.0513964336445958
    xp1, xp2 = -1.4210854715202004e-14, 20.000000000000007
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
  #=   saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
   #=  er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,60.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,60.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,60.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,60.0,1000.0)  
 =#
     ("B23",er1,er2)
   
end
err23=test23()
function test24()
    BSON.@load "sysB/formal_sysB24_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    u1, u2 = -0.03593893348744506, -1.739061066512555
    λ1, λ2 = -6.856244266050219, -0.04375573394978005
   
    c1, c2 = -2.6395715589141227, 0.6295715589141229
    xp1, xp2 = 0.0, 0.009999999999999998
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
   #=  saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2) =# 
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
  #=   er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,60.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,60.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,60.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,60.0,1000.0)   =#

     ("B24",er1,er2)
   
end
err24=test24()
function test25()
    BSON.@load "sysB/formal_sysB25_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    u1, u2 = -0.03593893348744506, -1.739061066512555
    λ1, λ2 = -6.856244266050219, -0.04375573394978005
   
    c1, c2 = -2.4761484322014207, -0.5238515677985796
    xp1, xp2 = -2.0, 1.0000000000000002
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
  #=   saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2) =# 
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
   #=  er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,60.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,60.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,60.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,60.0,1000.0)   =#

     ("B25",er1,er2)
   
end
err25=test25()
function test31()
    BSON.@load "sysB/formal_sysB31_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    u1, u2 = -0.02252283136759291, -1.6649771686324073
    λ1, λ2 = -6.909908674529628, -0.3400913254703717
   
    c1, c2 = -2.6362707559197682, 0.6362707559197684
    xp1, xp2 = 0.0, 0.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
   #=  saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
   #=  er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)  
 =#
     ("B31",er1,er2)
end
err31=test31()
function test32()
    BSON.@load "sysB/formal_sysB32_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    u1, u2 = -0.02252283136759291, -1.6649771686324073
    λ1, λ2 = -6.909908674529628, -0.3400913254703717
   
    c1, c2 = -3.6499836676620214, 0.6499836676620214
    xp1, xp2 = 1.1102230246251565e-16, 1.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
    #= saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
  #=   er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0) =#  

     ("B32",er1,er2)
   
end
err32=test32()
function test33()
    BSON.@load "sysB/formal_sysB33_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    u1, u2 = -0.02252283136759291, -1.6649771686324073
    λ1, λ2 = -6.909908674529628, -0.3400913254703717
   
    c1, c2 = -22.910528990764828, 0.910528990764828
    xp1, xp2 = 1.7763568394002505e-15, 20.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
   #=  saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
   #=  er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)   =#

     ("B33",er1,er2)
   
end
err33=test33()
function test34()
    BSON.@load "sysB/formal_sysB34_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    u1, u2 = -0.02252283136759291, -1.6649771686324073
    λ1, λ2 = -6.909908674529628, -0.3400913254703717
   
    c1, c2 = -2.646407885037191, 0.636407885037191
    xp1, xp2 = 8.673617379884035e-19, 0.01
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2

    #= saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
    #= er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)   =#

     ("B34",er1,er2)
   
end
err34=test34()
function test35()
    BSON.@load "sysB/formal_sysB35_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    u1, u2 = -0.02252283136759291, -1.6649771686324073
    λ1, λ2 = -6.909908674529628, -0.3400913254703717
   
    c1, c2 = -2.432293802791496, -0.5677061972085039
    xp1, xp2 = -2.0, 1.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2

    #= saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
   #=  er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)  =# 

     ("B35",er1,er2)
   
end
err35=test35()

function test41()
    BSON.@load "sysB/formal_sysB41_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    u1, u2 = -0.1144790629155213, -0.13539593708447864
    λ1, λ2 = -10.841674966758294, -9.168325033241706
   
    c1, c2 = -60.754387290570186, 58.754387290570186
    xp1, xp2 = 0.0, 0.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2

   #=  saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#

    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
   #=  er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)  
 =#
     ("B41",er1,er2)
   
end
err41=test41()
function test42()
    BSON.@load "sysB/formal_sysB42_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    u1, u2 = -0.1144790629155213, -0.13539593708447864
    λ1, λ2 = -10.841674966758294, -9.168325033241706
   
    c1, c2 = -67.22743560509412, 64.22743560509412
    xp1, xp2 = 0.0, 1.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2


    #= saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
   #=  er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)   =#

     ("B42",er1,er2)
   
end
err42=test42()

function test43()
    BSON.@load "sysB/formal_sysB43_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    u1, u2 = -0.1144790629155213, -0.13539593708447864
    λ1, λ2 = -10.841674966758294, -9.168325033241706
   
    c1, c2 = -190.21535358104902, 168.21535358104902
    xp1, xp2 = 0.0, 20.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2


  #=   saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
   #=  er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)   =#

     ("B43",er1,er2)
   
end
err43=test43()
 function test44()
    BSON.@load "sysB/formal_sysB44_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    u1, u2 = -0.1144790629155213, -0.13539593708447864
    λ1, λ2 = -10.841674966758294, -9.168325033241706
   
    c1, c2 = -60.81911777371541, 58.809117773715414
    xp1, xp2 = 0.0, 0.01
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2


  #=   saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
  #=   er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)  
 =#
     ("B44",er1,er2)
    
   
end
err44=test44()   =#

#= 
function test444()
    BSON.@load "sysB/formal_sysB44_mliqss2_errorE-3_step0.01_ft1000.bson" solmliqss2
    u1, u2 = -0.1144790629155213, -0.13539593708447864
    λ1, λ2 = -10.841674966758294, -9.168325033241706
   
    c1, c2 = -60.81911777371541, 58.809117773715414
    xp1, xp2 = 0.0, 0.01
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2


   #=  saveAbsoluteError(solmliqss2,1,x1)  
    saveAbsoluteError(solmliqss2,2,x2)  =#
    er1=getError(solmliqss2,1,x1)  
    er2=getError(solmliqss2,2,x2) 
  #=   er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)  
 =#
     ("B44",er1,er2)
   @show er1,er2
end
err44=test444()  =#
 function test45()
    BSON.@load "sysB/formal_sysB45_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    u1, u2 = -0.1144790629155213, -0.13539593708447864
    λ1, λ2 = -10.841674966758294, -9.168325033241706
   
    c1, c2 = 28.21275861611204, -31.18056545715832
    xp1, xp2 = -2.0, 1.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2


  #=   saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
   #=  er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)   =#

     ("B45",er1,er2)
   @show er1,er2
end
err45=test45() 


#= 


function test51()
    BSON.@load "sysB/formal_sysB51_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    u1, u2 = -8.73522174738572, -7.385745994549763
    λ1, λ2 = -10.841674966758294, -9.168325033241706
   
    c1, c2 = 11.68712513430148, -13.68712513430148
    xp1, xp2 = 0.0, 0.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2

    #= saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#

    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
   #=  er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)  
 =#
     ("B51",er1,er2)
   
   
end
err51=test51()
function test551()
    BSON.@load "sysB/formal_sysB51_mliqss2_errorE-3_step0.01_ft1000.bson" solmliqss2
    u1, u2 = -8.73522174738572, -7.385745994549763
    λ1, λ2 = -10.841674966758294, -9.168325033241706
   
    c1, c2 = 11.68712513430148, -13.68712513430148
    xp1, xp2 = 0.0, 0.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2

  #=   saveAbsoluteError(solmliqss2,1,x1)  
    saveAbsoluteError(solmliqss2,2,x2)  =#

    er1=getError(solmliqss2,1,x1)  
    er2=getError(solmliqss2,2,x2) 
   #=  er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)  
 =#
     ("B51",er1,er2)
    
   
end
err51=test551()
 function test52()
    BSON.@load "sysB/formal_sysB52_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    u1, u2 = -8.73522174738572, -7.385745994549763
    λ1, λ2 = -10.841674966758294, -9.168325033241706
   
    c1, c2 = 17.16017344882542, -20.16017344882542
    xp1, xp2 = 0.0, 1.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2


  #=   saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
   #=  er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)   =#

     ("B52",er1,er2)
   
end
err52=test52()

function test53()
    BSON.@load "sysB/formal_sysB53_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    u1, u2 = -8.73522174738572, -7.385745994549763
    λ1, λ2 = -10.841674966758294, -9.168325033241706
   
    c1, c2 = 121.14809142478035, -143.14809142478035
    xp1, xp2 = 0.0, 20.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2


    #= saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2) =# 
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
   #=  er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)   =#

     ("B53",er1,er2)
   
end
err53=test53()
function test54()
    BSON.@load "sysB/formal_sysB54_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    u1, u2 = -8.73522174738572, -7.385745994549763
    λ1, λ2 = -10.841674966758294, -9.168325033241706
   
    c1, c2 = 11.74185561744672, -13.75185561744672
    xp1, xp2 = 0.0, 0.01
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2


   #=  saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
  #=   er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)  
 =#
     ("B54",er1,er2)
   
end
err54=test54() 
function test55()
    BSON.@load "sysB/formal_sysB55_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    u1, u2 = -8.73522174738572, -7.385745994549763
    λ1, λ2 = -10.841674966758294, -9.168325033241706
   
    c1, c2 = 15.67811643831823, -18.67811643831823
    xp1, xp2 = -2.0, 1.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2


    #= saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
   #=  er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)   =#

     ("B55",er1,er2)
   
end
err55=test55()

 =#





#= 




 XLSX.openxlsx("formal_Interp_______.xlsx", mode="w") do xf
    sheet = xf[1]
   
    sheet["A1"] = "formal_interp"
   
    sheet["A4"] = collect(("Sys","x1","x2"))
    sheet["A5"] = collect(err11)
    sheet["A6"] = collect(err12)
    sheet["A7"] = collect(err13)
    sheet["A8"] = collect(err14)
    sheet["A9"] = collect(err15)

    sheet["A10"] = collect(err21)
    sheet["A11"] = collect(err22)
    sheet["A12"] = collect(err23)
    sheet["A13"] = collect(err24)
    sheet["A14"] = collect(err25)

    sheet["A15"] = collect(err31)
    sheet["A16"] = collect(err32)
    sheet["A17"] = collect(err33)
    sheet["A18"] = collect(err34)
    sheet["A19"] = collect(err35)

    sheet["A20"] = collect(err41)
    sheet["A21"] = collect(err42)
    sheet["A22"] = collect(err43)
    sheet["A23"] = collect(err44) 
    sheet["A24"] = collect(err45)
    sheet["A25"] = collect(err51)
    sheet["A26"] = collect(err52)
    sheet["A27"] = collect(err53)
    sheet["A28"] = collect(err54)
    sheet["A29"] = collect(err55)

   end  =#