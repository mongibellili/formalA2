
using XLSX
using relaxedqssA
using StaticArrays
using BSON

#include("/home/unknown/formalqssAB/formalqssAB/src/models/InvestigateBProblem.jl") 
 function test11()
    BSON.@load "sysA/relaxedqssA_sys1_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
   
   
    u1, u2 = -0.5000000000000002 + 0.8660254037844386im, -0.5000000000000002 - 0.8660254037844386im
    λ1, λ2 = -0.5 - 0.8660254037844385im, -0.5 + 0.8660254037844385im
   
    c1, c2 = 0.0 - 0.5773502691896258im, 0.0 + 0.5773502691896258im
    xp1, xp2 = 0.0, 0.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
  #  saveAbsoluteError(solmliqss2Interp,1,x1)  
   # saveAbsoluteError(solmliqss2Interp,2,x2) 
    er1=Float64(getError(solmliqss2Interp,1,x1))  
    er2=Float64(getError(solmliqss2Interp,2,x2)) 
   #=  er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)  =# 



     ("A1",er1,er2)
   
end
err11=test11()
function test12()
    BSON.@load "sysA/relaxedqssA_sys2_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    u1, u2 = -0.00010001000200050016, -0.9998999899979997
    λ1, λ2 = -99.98999899979995, -0.010001000200048793
   
    c1, c2 = -19.203841152384133, 19.203841152384133
    xp1, xp2 = 20.2, 0.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
   # saveAbsoluteError(solmliqss2Interp,1,x1)  
   # saveAbsoluteError(solmliqss2Interp,2,x2) 
    er1=Float64(getError(solmliqss2Interp,1,x1))  
    er2=Float64(getError(solmliqss2Interp,2,x2)) 
    #= er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)   =#

     ("A2",er1,er2)
   
end
err12=test12()
function test13()
    BSON.@load "sysA/relaxedqssA_sys3_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    u1, u2 = 0.0 - 1.0im, -0.0 + 1.0im
    λ1, λ2 = -1.0 - 1.0im, -1.0 + 1.0im
   
    c1, c2 = 1.65 - 1.75im, 1.65 + 1.75im
    xp1, xp2 = -0.5, 0.7
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
   # saveAbsoluteError(solmliqss2Interp,1,x1)  
   # saveAbsoluteError(solmliqss2Interp,2,x2) 
    er1=Float64(getError(solmliqss2Interp,1,x1))  
    er2=Float64(getError(solmliqss2Interp,2,x2)) 
    #= er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)   =#

     ("A3",er1,er2)
   
end
err13=test13()
function test14()
    BSON.@load "sysA/relaxedqssA_sys4_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    u1, u2 = -0.0 - 0.9999999999999998im, -0.0 + 0.9999999999999998im
    λ1, λ2 = -1.0 - 10.000000000000002im, -1.0 + 10.000000000000002im
   
    c1, c2 = 1.9841584158415841 - 1.941584158415842im, 1.9841584158415841 + 1.941584158415842im
    xp1, xp2 = -0.11683168316831682, 0.03168316831683169
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
   # saveAbsoluteError(solmliqss2Interp,1,x1)  
   # saveAbsoluteError(solmliqss2Interp,2,x2) 
    er1=Float64(getError(solmliqss2Interp,1,x1))  
    er2=Float64(getError(solmliqss2Interp,2,x2)) 
  #=   er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)  =# 

     ("A4",er1,er2)
end
err14=test14()
function test15()
    BSON.@load "sysA/relaxedqssA_sys5_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    u1, u2 = -3.8789024427351744, 0.12890244273517468
    λ1, λ2 = -4.128902442735175, -0.12109755726482532
   
    c1, c2 = -24.89283101895126, 20.89283101895126
    xp1, xp2 = 0.75, 4.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
   # saveAbsoluteError(solmliqss2Interp,1,x1)  
   # saveAbsoluteError(solmliqss2Interp,2,x2) 
    er1=Float64(getError(solmliqss2Interp,1,x1))  
    er2=Float64(getError(solmliqss2Interp,2,x2)) 
   #=  er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)   =#

     ("A5",er1,er2)
   
end
err15=test15()

function test16()
    BSON.@load "sysA/relaxedqssA_sys6_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    u1, u2 = -0.0004900990099009912 - 3.1465838394698573im, -0.0004900990099009912 + 3.1465838394698573im
    λ1, λ2 = -1.0050499999999998 - 31.780496778645546im, -1.0050499999999998 + 31.780496778645546im
   
    c1, c2 = -3.846103396992669 + 0.004227711787777913im, -3.846103396992669 - 0.004227711787777913im
    xp1, xp2 = 0.9696243578875977, 9.992206793985337
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
   # saveAbsoluteError(solmliqss2Interp,1,x1)  
   # saveAbsoluteError(solmliqss2Interp,2,x2) 
    er1=Float64(getError(solmliqss2Interp,1,x1))  
    er2=Float64(getError(solmliqss2Interp,2,x2)) 
    #= er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,60.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,60.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,60.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,60.0,1000.0)   =#

     ("A6",er1,er2)
   
end
err16=test16()


 XLSX.openxlsx("relaxedqssA_Interp_sysA.xlsx", mode="w") do xf
    sheet = xf[1]
   
    sheet["A1"] = "relaxedqssA_interp"
   
    sheet["A4"] = collect(("SysA","x1","x2"))
     sheet["A5"] = collect(err11)
    sheet["A6"] = collect(err12)
    sheet["A7"] = collect(err13)
    sheet["A8"] = collect(err14)
    sheet["A9"] = collect(err15)

    sheet["A10"] = collect(err16)
   

   end 