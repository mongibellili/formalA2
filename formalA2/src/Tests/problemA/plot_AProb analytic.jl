#using qssmodifiedv01
using Plots
#include("D://QS_Solver//qss//src//models//InvestigateBProblem.jl") 
 function test11()
  
   
   
    u1, u2 = -0.5000000000000002 + 0.8660254037844386im, -0.5000000000000002 - 0.8660254037844386im
    λ1, λ2 = -0.5 - 0.8660254037844385im, -0.5 + 0.8660254037844385im
   
    c1, c2 = 0.0 - 0.5773502691896258im, 0.0 + 0.5773502691896258im
    xp1, xp2 = 0.0, 0.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
    p1=plot()
    p1=plot!(p1,x1,title="A1",label="x1",xlims=(0,100),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_A1 Analytic_FT100.png")
   
     
   
end
test11()
function test12()
   
    u1, u2 = -0.00010001000200050016, -0.9998999899979997
    λ1, λ2 = -99.98999899979995, -0.010001000200048793
   
    c1, c2 = -19.203841152384133, 19.203841152384133
    xp1, xp2 = 20.2, 0.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
    
    p1=plot()
    p1=plot!(p1,x1,title="A2",label="x1",xlims=(0,100),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_A2 Analytic_FT100.png")
end
test12()
function test13()
  
    u1, u2 = 0.0 - 1.0im, -0.0 + 1.0im
    λ1, λ2 = -1.0 - 1.0im, -1.0 + 1.0im
   
    c1, c2 = 1.65 - 1.75im, 1.65 + 1.75im
    xp1, xp2 = -0.5, 0.7
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
    
    p1=plot()
    p1=plot!(p1,x1,title="A3",label="x1",xlims=(0,100),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_A3 Analytic_FT100.png")
   
end
test13()
function test14()
   
    u1, u2 = -0.0 - 0.9999999999999998im, -0.0 + 0.9999999999999998im
    λ1, λ2 = -1.0 - 10.000000000000002im, -1.0 + 10.000000000000002im
   
    c1, c2 = 1.9841584158415841 - 1.941584158415842im, 1.9841584158415841 + 1.941584158415842im
    xp1, xp2 = -0.11683168316831682, 0.03168316831683169
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
    
    p1=plot()
    p1=plot!(p1,x1,title="A4",label="x1",xlims=(0,100),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_A4 Analytic_FT100.png")
   
end
test14()
function test15()
   
    u1, u2 = -3.8789024427351744, 0.12890244273517468
    λ1, λ2 = -4.128902442735175, -0.12109755726482532
   
    c1, c2 = -24.89283101895126, 20.89283101895126
    xp1, xp2 = 0.75, 4.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
    
    p1=plot()
    p1=plot!(p1,x1,title="A5",label="x1",xlims=(0,100),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_A5 Analytic_FT100.png")
   
end
test15()
function test16()
   
    u1, u2 = -0.0004900990099009912 - 3.1465838394698573im, -0.0004900990099009912 + 3.1465838394698573im
    λ1, λ2 = -1.0050499999999998 - 31.780496778645546im, -1.0050499999999998 + 31.780496778645546im
   
    c1, c2 = -3.846103396992669 + 0.004227711787777913im, -3.846103396992669 - 0.004227711787777913im
    xp1, xp2 = 0.9696243578875977, 9.992206793985337
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
    
    p1=plot()
    p1=plot!(p1,x1,title="A6",label="x1",xlims=(0,100),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_A6 Analytic_FT100.png")
   
end
test16()
