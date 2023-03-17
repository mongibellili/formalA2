
using XLSX
using relaxedqssA
using StaticArrays
using BSON

#include("/home/unknown/formalqssAB/formalqssAB/src/models/InvestigateBProblem.jl") 
function test11()
    BSON.@load "sys1/relaxedqssA_sys12_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
   
   
    u1=(7*sqrt(51)-50)/100
   u2=(-7*sqrt(51)-50)/100
   λ1=-7*sqrt(51)-50
   λ2=7*sqrt(51)-50
   c2=960/(7*sqrt(51))
   c1=-c2
   
   x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+20.2
   x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t) 
    saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2) 
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
   #=  er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)  =# 

     ("B11",er1,er2)
   
end
err11=test11()







 XLSX.openxlsx("relaxedqssA_sys12_Interp.xlsx", mode="w") do xf
    sheet = xf[1]
   
    sheet["A1"] = "relaxedqssA_sys12_interp"
   
    sheet["A4"] = collect(("Sys","x1","x2"))
    sheet["A5"] = collect(err11)
    

   end