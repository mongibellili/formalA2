#= #= using Symbolics

 @variables t x y
 D = Differential(t)

z = t + t^2
@show  D(z) # symbolic representation of derivative(t + t^2, t)
 =#
using Plots
using Polynomials
#= bi=1.0;ci=1.1;βi=2.03;αi=-4.015
 f(y)=(bi*y+ci*y*y)/(1+αi*y+βi*y*y)-5.0
# p = Polynomial([βi, αi, 1])
 #@show roots(p)
 p1=plot!(f,xlim=[0.0,3.0],ylim=[-10.0,10.0]) =#

 #= bi=1.0;ci=1.1;βi=22.03;αi=-15.015
 g(y)=5.0+(bi*y+ci*y*y)/(1+αi*y+βi*y*y)
 #=  p = Polynomial([1, αi, βi])
 @show roots(p) =#

 p1=plot!(g,xlim=[0.0,3.0],ylim=[-20.0,20.0])


 savefig(p1, "plot_fj.png") =#


 #bi=0.004;ci=-0.005;βi=7.75;αi=6.0
 #= bi=1.0;ci=1.1;βi=2.03;αi=-4.015
 g(y)=(bi*y+ci*y*y)/(1+αi*y+βi*y*y)
  p = Polynomial([1, αi, βi])
 #@show roots(p)

 p1=plot!(g,xlim=[-3.0,6.0],ylim=[-10.0,10.0])


 savefig(p1, "plot_fi.png") =#

 #bi=0.90;ci=-0.6;βi=2.03;αi=-4.015
 αi, βi, quani, ci, bi, resi1, resi2, resi3, resi4 = 25.5, 32.5, 1.0e-6, -0.007060038244102174, 0.005844076254252428, 0.8202080625013038, 0.00017189936790235752, 0.8353954241513185, -0.00017033531193132506
 g(y)=(bi*y+ci*y*y)/(1+αi*y+βi*y*y)
  #= p = Polynomial([1, αi, βi])
 @show roots(p) =#

 p1=plot!(g,xlim=[-0.01,30.0],ylim=[-1.0e-2,1.5e-5])
 delta(y)=1.0e-6
 p1=plot!(delta,xlim=[-0.01,30.0],ylim=[-1.0e-2,1.5e-5])
 delta_(y)=-1.0e-6
 p1=plot!(delta_,)
 savefig(p1, "plot_fl6.png")

 #= bi=-0.004;ci=-0.005;βi=7.75;αi=6.0
 g(y)=(bi*y+ci*y*y)


 p1=plot!(g,xlim=[-1.0,0.0],ylim=[0.0005,1e-3])


 savefig(p1, "plot_f2.png") =#
 =#
 
 using LinearAlgebra
 A=[-0.75 0.2;-12.65 2.0]
#=  V=eigvecs(A)
 λ=eigvals(A)
 #λ = ComplexF64[0.625 - 0.7996092795859737im, 0.625 + 0.7996092795859737im]
#V = ComplexF64[0.10784645438255967 + 0.06271638232344737im 0.10784645438255967 - 0.06271638232344737im; 0.992187380319549 - 0.0im 0.992187380319549 + 0.0im]
 @show λ
 @show V =#

 Δ_=sqrt(-((A[1,1]-A[2,2])^2+4*A[1,2]*A[2,1]))
 α=-A[1,1]-A[2,2]
 λ=(-α-Δ_*im)/2
 #@show λ 
 n1=(-A[1,1]+A[2,2]+im*Δ_)/(-2*A[2,1])
 @show n1

 
 #=  f=(0.21+sqrt(1.2589))/(2)
    @show f =#