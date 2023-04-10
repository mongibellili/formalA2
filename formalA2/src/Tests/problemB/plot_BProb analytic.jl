#using qssmodifiedv01
using Plots
#include("D://QS_Solver//qss//src//models//InvestigateBProblem.jl") 
#= function test11()
  
   
   
    u1, u2 = -0.06805253887502617, -0.5248046039821167
    λ1, λ2 = -3.523632227874817, -0.32636777212518275
   
    c1, c2 = -4.487356192869495, 2.4873561928694947
    xp1, xp2 = 0.0, 0.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
    p1=plot()
    p1=plot!(p1,x1,title="B11",label="x1",xlims=(0,1000),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_B11 Analytic_FT1000.png")
   
     
   
end
test11()
function test12()
   
    u1, u2 = -0.06805253887502617, -0.5248046039821167
    λ1, λ2 = -3.523632227874817, -0.32636777212518275
   
    c1, c2 = -5.636348488851935, 2.6363484888519344
    xp1, xp2 = 1.1102230246251565e-16, 0.9999999999999999
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
    
    p1=plot()
    p1=plot!(p1,x1,title="B12",label="x1",xlims=(0,1000),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_B12 Analytic_FT1000.png")
end
test12()
function test13()
  
    u1, u2 = -0.06805253887502617, -0.5248046039821167
    λ1, λ2 = -3.523632227874817, -0.32636777212518275
   
    c1, c2 = -27.467202112518294, 5.467202112518295
    xp1, xp2 = 3.552713678800501e-15, 20.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
    
    p1=plot()
    p1=plot!(p1,x1,title="B13",label="x1",xlims=(0,1000),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_B13 Analytic_FT1000.png")
   
end
test13()
function test14()
   
    u1, u2 = -0.06805253887502617, -0.5248046039821167
    λ1, λ2 = -3.523632227874817, -0.32636777212518275
   
    c1, c2 = -4.498846115829319, 2.488846115829319
    xp1, xp2 = 0.0, 0.009999999999999998
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
    
    p1=plot()
    p1=plot!(p1,x1,title="B14",label="x1",xlims=(0,1000),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_B14 Analytic_FT1000.png")
   
end
test14()
function test15()
   
    u1, u2 = -0.06805253887502617, -0.5248046039821167
    λ1, λ2 = -3.523632227874817, -0.32636777212518275
   
    c1, c2 = -1.2576052870427024, -1.7423947129572968
    xp1, xp2 = -2.0, 0.9999999999999993
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
    
    p1=plot()
    p1=plot!(p1,x1,title="B15",label="x1",xlims=(0,1000),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_B15 Analytic_FT1000.png")
   
end
test15()
function test21()
   
    u1, u2 = -0.03593893348744506, -1.739061066512555
    λ1, λ2 = -6.856244266050219, -0.04375573394978005
   
    c1, c2 = -2.6293605409677845, 0.6293605409677845
    xp1, xp2 = 0.0, 0.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
    
    p1=plot()
    p1=plot!(p1,x1,title="B21",label="x1",xlims=(0,1000),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_B21 Analytic_FT1000.png")
   
end
test21()
function test22()
    
    u1, u2 = -0.03593893348744506, -1.739061066512555
    λ1, λ2 = -6.856244266050219, -0.04375573394978005
   
    c1, c2 = -3.6504623356016257, 0.6504623356016255
    xp1, xp2 = 0.0, 1.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
    
    p1=plot()
    p1=plot!(p1,x1,title="B22",label="x1",xlims=(0,1000),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_B22 Analytic_FT1000.png")
end
test22()
function test23()
    
    u1, u2 = -0.03593893348744506, -1.739061066512555
    λ1, λ2 = -6.856244266050219, -0.04375573394978005
   
    c1, c2 = -23.0513964336446, 1.0513964336445958
    xp1, xp2 = -1.4210854715202004e-14, 20.000000000000007
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
    
    p1=plot()
    p1=plot!(p1,x1,title="B23",label="x1",xlims=(0,1000),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_B23 Analytic_FT1000.png")
end
test23()
function test24()
   
    u1, u2 = -0.03593893348744506, -1.739061066512555
    λ1, λ2 = -6.856244266050219, -0.04375573394978005
   
    c1, c2 = -2.6395715589141227, 0.6295715589141229
    xp1, xp2 = 0.0, 0.009999999999999998
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
    p1=plot()
    p1=plot!(p1,x1,title="B24",label="x1",xlims=(0,1000),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_B24 Analytic_FT1000.png")
end
test24()
function test25()
   
    u1, u2 = -0.03593893348744506, -1.739061066512555
    λ1, λ2 = -6.856244266050219, -0.04375573394978005
   
    c1, c2 = -2.4761484322014207, -0.5238515677985796
    xp1, xp2 = -2.0, 1.0000000000000002
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
    p1=plot()
    p1=plot!(p1,x1,title="B25",label="x1",xlims=(0,1000),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_B25 Analytic_FT1000.png") 
end
test25()
function test31()
   
    u1, u2 = -0.02252283136759291, -1.6649771686324073
    λ1, λ2 = -6.909908674529628, -0.3400913254703717
   
    c1, c2 = -2.6362707559197682, 0.6362707559197684
    xp1, xp2 = 0.0, 0.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
    p1=plot()
    p1=plot!(p1,x1,title="B31",label="x1",xlims=(0,1000),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_B31 Analytic_FT1000.png")
end
test31()
function test32()
   
    u1, u2 = -0.02252283136759291, -1.6649771686324073
    λ1, λ2 = -6.909908674529628, -0.3400913254703717
   
    c1, c2 = -3.6499836676620214, 0.6499836676620214
    xp1, xp2 = 1.1102230246251565e-16, 1.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
    p1=plot()
    p1=plot!(p1,x1,title="B32",label="x1",xlims=(0,1000),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_B32 Analytic_FT1000.png")
end
test32()
function test33()
    
    u1, u2 = -0.02252283136759291, -1.6649771686324073
    λ1, λ2 = -6.909908674529628, -0.3400913254703717
   
    c1, c2 = -22.910528990764828, 0.910528990764828
    xp1, xp2 = 1.7763568394002505e-15, 20.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
    p1=plot()
    p1=plot!(p1,x1,title="B33",label="x1",xlims=(0,1000),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_B33 Analytic_FT1000.png")
end
test33()
function test34()
   
    u1, u2 = -0.02252283136759291, -1.6649771686324073
    λ1, λ2 = -6.909908674529628, -0.3400913254703717
   
    c1, c2 = -2.646407885037191, 0.636407885037191
    xp1, xp2 = 8.673617379884035e-19, 0.01
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
    p1=plot()
    p1=plot!(p1,x1,title="B34",label="x1",xlims=(0,1000),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_B34 Analytic_FT1000.png")
end
test34()
function test35()
    
    u1, u2 = -0.02252283136759291, -1.6649771686324073
    λ1, λ2 = -6.909908674529628, -0.3400913254703717
   
    c1, c2 = -2.432293802791496, -0.5677061972085039
    xp1, xp2 = -2.0, 1.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
    p1=plot()
    p1=plot!(p1,x1,title="B35",label="x1",xlims=(0,1000),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_B35 Analytic_FT1000.png")
end
test35() 
=#
#= function test41()
    u1, u2 = -0.1144790629155213, -0.13539593708447864
    λ1, λ2 = -10.841674966758294, -9.168325033241706
    c1, c2 = -60.754387290570186, 58.754387290570186
    xp1, xp2 = 0.0, 0.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
    p1=plot()
    p1=plot!(p1,x1,title="B41",label="x1",xlims=(0,100),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_B41 Analytic_FT100.png")
    #= @show x2(0.01)
    @show x2(0.1)
    @show x2(1.0) =#
end
test41()

 function test42()
   
    u1, u2 = -0.1144790629155213, -0.13539593708447864
    λ1, λ2 = -10.841674966758294, -9.168325033241706
   
    c1, c2 = -67.22743560509412, 64.22743560509412
    xp1, xp2 = 0.0, 1.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
    
    p1=plot()
    p1=plot!(p1,x1,title="B42",label="x1",xlims=(0,100),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_B42 Analytic_FT100.png")
end
test42() 

function test43()
  
    u1, u2 = -0.1144790629155213, -0.13539593708447864
    λ1, λ2 = -10.841674966758294, -9.168325033241706
   
    c1, c2 = -190.21535358104902, 168.21535358104902
    xp1, xp2 = 0.0, 20.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
    
    p1=plot()
    p1=plot!(p1,x1,title="B43",label="x1",xlims=(0,100),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_B43 Analytic_FT100.png")
   
end
test43()=#
#= function test44()
    u1, u2 = -0.1144790629155213, -0.13539593708447864
    λ1, λ2 = -10.841674966758294, -9.168325033241706
   
    c1, c2 = -60.81911777371541, 58.809117773715414
    xp1, xp2 = 0.0, 0.01
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
  
    p1=plot()
    p1=plot!(p1,x1,title="B44",label="x1",xlims=(0,124.9),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_B44 Analytic_FT124_9.png")
    i=0.0
    while i<1
        println(x2(i))
        i=i+0.01
    end
end
test44()  =# 
#=  function test45()
   
    u1, u2 = -0.1144790629155213, -0.13539593708447864
    λ1, λ2 = -10.841674966758294, -9.168325033241706
   
    c1, c2 = 28.21275861611204, -31.18056545715832
    xp1, xp2 = -2.0, 1.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
    
    p1=plot()
    p1=plot!(p1,x1,title="B45",label="x1",xlims=(0,100),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_B45 Analytic_FT100.png")
   
end
test45()  =#


 function test51()
    u1, u2 = -8.73522174738572, -7.385745994549763
    λ1, λ2 = -10.841674966758294, -9.168325033241706
    c1, c2 = 11.68712513430148, -13.68712513430148
    xp1, xp2 = 0.0, 0.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
    p1=plot()
    p1=plot!(p1,x1,title="B51",label="x1",xlims=(0,0.5),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_B51 Analytic_FT100.png")
   #=  @show x2(0.01)
    @show x2(0.1)
    @show x2(1.0) =#
end
test51()
#=
function test52()
   
    u1, u2 = -8.73522174738572, -7.385745994549763
    λ1, λ2 = -10.841674966758294, -9.168325033241706
   
    c1, c2 = 17.16017344882542, -20.16017344882542
    xp1, xp2 = 0.0, 1.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
    
    p1=plot()
    p1=plot!(p1,x1,title="B52",label="x1",xlims=(0,100),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_B52 Analytic_FT100.png")
end
test52() 

function test53()
  
    u1, u2 = -8.73522174738572, -7.385745994549763
    λ1, λ2 = -10.841674966758294, -9.168325033241706
   
    c1, c2 = 121.14809142478035, -143.14809142478035
    xp1, xp2 = 0.0, 20.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
    
    p1=plot()
    p1=plot!(p1,x1,title="B53",label="x1",xlims=(0,100),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_B53 Analytic_FT100.png")
   
end
test53()
function test54()
    u1, u2 = -8.73522174738572, -7.385745994549763
    λ1, λ2 = -10.841674966758294, -9.168325033241706
   
    c1, c2 = 11.74185561744672, -13.75185561744672
    xp1, xp2 = 0.0, 0.01
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
    
    p1=plot()
    p1=plot!(p1,x1,title="B54",label="x1",xlims=(0,100),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_B54 Analytic_FT100.png")
   
end
test54() 
 function test55()
   
    u1, u2 = -8.73522174738572, -7.385745994549763
    λ1, λ2 = -10.841674966758294, -9.168325033241706
   
    c1, c2 = 15.67811643831823, -18.67811643831823
    xp1, xp2 = -2.0, 1.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
    
    p1=plot()
    p1=plot!(p1,x1,title="B55",label="x1",xlims=(0,100),)
    p1=plot!(p1,x2,label="x2")
    savefig(p1, "plot_B55 Analytic_FT100.png")
   
end
test55()=#
 