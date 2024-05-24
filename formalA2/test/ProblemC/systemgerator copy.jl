using LinearAlgebra
struct LTIProblem
    name::String
    x0::Vector{Float64}#initconds
    input::Vector{Float64}
    jac::Array{Float64, 2}
end
struct LTISolution
    name::String
    eignVal::Vector{Float64}
    eignVec::Vector{Float64}
    solPart::Vector{Float64}
    coefs::Vector{Float64}
end
struct LTIComplexSolution
    name::String
    solPart::Vector{Float64}
    coefs::Vector{Float64}
    delta::Float64
end

function generateSol(ltiprob::LTIProblem)
    A=ltiprob.jac;b=ltiprob.input;x10=ltiprob.x0[1];x20=ltiprob.x0[2]
    xp=-inv(A)*b 
    V=eigvecs(A)
    λ=eigvals(A)
   # @show λ[1],λ[2]
    V1=V[1]/V[2] 
    V2=V[3]/V[4] 
  #  @show V1,V2
    c2=(x10-xp[1]-(x20-xp[2])*V1)/(V2-V1)
    c1=x20-xp[2]-c2
   #=  @show xp[1],xp[2]
    @show c1,c2 =#
    LTISolution(ltiprob.name,[λ[1],λ[2]],[V1,V2],[xp[1],xp[2]],[c1,c2])
end
function generateSolComplex(ltiprob::LTIProblem)
    A=ltiprob.jac;b=ltiprob.input;x10=ltiprob.x0[1];x20=ltiprob.x0[2]
    xp=-inv(A)*b 
  #@show xp
  #=   V=eigvecs(A)
    λ=eigvals(A)
   # @show λ[1],λ[2]
    V1=V[1]/V[2] 
    V2=V[3]/V[4]  =#
  #  @show V1,V2
  Δ_=sqrt(-((A[1,1]-A[2,2])^2+4*A[1,2]*A[2,1]))
  α=-A[1,1]-A[2,2]
  c1=x20-xp[2]
 
  c2=((-2*A[2,1])*(x10-xp[1])-(x20-xp[2])*(A[2,2]-A[1,1]))/Δ_
    
   #=  @show xp[1],xp[2]
    @show c1,c2 =#
    LTIComplexSolution(ltiprob.name,[xp[1],xp[2]],[c1,c2],Δ_)
end

function consructJacobiansSystemAA(coeffs:: Vector{Float64},posCoeffs:: Vector{Float64},negCoeffs:: Vector{Float64})
    allJacsys=Vector{ Array{Float64, 2}}()
    for a11 in coeffs
           for a22 in coeffs
               if a11+a22<0
                   for a12 in coeffs
                       for a21 in coeffs
                             if a11*a22!=a12*a21 # to have A invertible
                               if (a11-a22)^2+4*(a12*a21) >0 # this is just to have real solutions
                                   A=[a11 a12;a21 a22]
                                   push!(allJacsys,A)
                               end
                            end
                         
                       end
                   end
               end
           end
     
   end
  #  push!(allJacsys,A)
   # display(allJacsys)
   return allJacsys
end
function consructJacobiansSystemBB(coeffs:: Vector{Float64},posCoeffs:: Vector{Float64},negCoeffs:: Vector{Float64})
    allJacsys=Vector{ Array{Float64, 2}}()
    for a11 in coeffs
           for a22 in coeffs
               if a11+a22>0
                   for a12 in coeffs
                       for a21 in coeffs
                             if a11*a22!=a12*a21 # to have A invertible
                               if (a11-a22)^2+4*(a12*a21) >0 # this is just to have real solutions
                                   A=[a11 a12;a21 a22]
                                   push!(allJacsys,A)
                               end
                            end
                         
                       end
                   end
               end
           end
     
   end
  #  push!(allJacsys,A)
   # display(allJacsys)
   return allJacsys
end

function consructJacobiansSystemCC(coeffs:: Vector{Float64},posCoeffs:: Vector{Float64},negCoeffs:: Vector{Float64})
    allJacsys=Vector{ Array{Float64, 2}}()
    for a11 in coeffs
           for a22 in coeffs
               if a11+a22<0
                   for a12 in coeffs
                       for a21 in coeffs
                             if a11*a22!=a12*a21 # to have A invertible
                               if (a11-a22)^2+4*(a12*a21) <0 # this is just to have real solutions
                                   A=[a11 a12;a21 a22]
                                   push!(allJacsys,A)
                               end
                            end
                         
                       end
                   end
               end
           end
     
   end
  #  push!(allJacsys,A)
   # display(allJacsys)
   return allJacsys
end
function consructJacobiansSystemDD(coeffs:: Vector{Float64},posCoeffs:: Vector{Float64},negCoeffs:: Vector{Float64})
    allJacsys=Vector{ Array{Float64, 2}}()
    for a11 in coeffs
           for a22 in coeffs
               if a11+a22>0
                   for a12 in coeffs
                       for a21 in coeffs
                             if a11*a22!=a12*a21 # to have A invertible
                               if (a11-a22)^2+4*(a12*a21) <0 # this is just to have real solutions
                                   A=[a11 a12;a21 a22]
                                   push!(allJacsys,A)
                               end
                            end
                         
                       end
                   end
               end
           end
     
   end
  #  push!(allJacsys,A)
   # display(allJacsys)
   return allJacsys
end

function writeProblemstoFile(prbName,coefs,posCoefs,negCoefs,sysFunc,initConds,us,ul,counter)
    jacs=sysFunc(coefs,posCoefs,negCoefs)
    inputs=[[us,us],[-us,-us],[-us,us],[us,ul],[-us,-ul],[-us,ul],[us,-ul],[ul,ul],[-ul,-ul],[-ul,ul]]
    #inputs=[[us,us],[us,-ul]]
   #@show jacs
    #counter=0
 
    path="./type$(prbName).jl" #default path
    #@show path
    vectprStr="["  # helper to manually written vector later in maintest
    for jac in jacs
        for input in inputs
            counter+=1
            lti_prob=LTIProblem("$(prbName)_$(counter)",initConds,input,jac)
            anaSol=generateSol(lti_prob)
           # push!(allCsols,anaSol)
           vectprStr*="$(prbName)_$(counter),"
            ss="\n function $(lti_prob.name)() \n"
            ss*=" odeprob = @NLodeProblem begin
                name=($(lti_prob.name),)
                u = [$(lti_prob.x0[1]), $(lti_prob.x0[2])]
                du[1] = $(lti_prob.jac[1,1])*u[1]+$(lti_prob.jac[1,2])*u[2]+$(lti_prob.input[1])
                du[2] =$(lti_prob.jac[2,1])*u[1]+$(lti_prob.jac[2,2])*u[2]+$(lti_prob.input[2])
            end  "
            ss*="\n x1(t)=$(anaSol.coefs[1])*$(anaSol.eignVec[1])*exp($(anaSol.eignVal[1])*t)+$(anaSol.coefs[2])*$(anaSol.eignVec[2])*exp($(anaSol.eignVal[2])*t)+$(anaSol.solPart[1])"
            ss*="\n x2(t)=$(anaSol.coefs[1])*exp($(anaSol.eignVal[1])*t)+$(anaSol.coefs[2])*exp($(anaSol.eignVal[2])*t)+$(anaSol.solPart[2])"
            ss*="\n return (odeprob,x1,x2) \n end "
            open(path, "a") do io    
                println(io,ss) 
            end
        end
    end
    vectprStr*="]"
    open(path, "a") do io    
        println(io,vectprStr) 
    end
end


function writeComplexProblemstoFile(prbName,coefs,posCoefs,negCoefs,sysFunc,initConds,us,ul,counter)
    jacs=sysFunc(coefs,posCoefs,negCoefs)
    #jacs=[[3.0 -13.0;5.0 1.0]]
    #inputs=[[us,us],[-us,-us],[-us,us],[us,ul],[-us,-ul],[-us,ul],[us,-ul],[ul,ul],[-ul,-ul],[-ul,ul]]
    inputs=[[us,us],[-us,ul]]
   #@show jacs
    #counter=0
 
    path="./type$(prbName).jl" #default path
   # @show path
    vectprStr="["  # helper to manually written vector later in maintest
    for jac in jacs
        for input in inputs
            counter+=1
            lti_prob=LTIProblem("$(prbName)_$(counter)",initConds,input,jac)
            anaSol=generateSolComplex(lti_prob)
           # push!(allCsols,anaSol)
           sqΔ=anaSol.delta
           α=-jac[1,1]-jac[2,2]
           aiijj=jac[2,2]-jac[1,1]
           vectprStr*="$(prbName)_$(counter),"
            ss="\n function $(lti_prob.name)() \n"
            ss*=" odeprob = @NLodeProblem begin
                name=($(lti_prob.name),)
                u = [$(lti_prob.x0[1]), $(lti_prob.x0[2])]
                du[1] = $(lti_prob.jac[1,1])*u[1]+$(lti_prob.jac[1,2])*u[2]+$(lti_prob.input[1])
                du[2] =$(lti_prob.jac[2,1])*u[1]+$(lti_prob.jac[2,2])*u[2]+$(lti_prob.input[2])
            end  "
            ss*="\n x1(t)=$(-anaSol.coefs[1]/(2*jac[2,1]))*exp($(-α/2)*t)*($(aiijj)*cos($(sqΔ/2)*t)+$(sqΔ)*sin($(sqΔ/2)*t))+$(-anaSol.coefs[2]/(2*jac[2,1]))*exp($(-α/2)*t)*($(sqΔ)*cos($(sqΔ/2)*t)+$(-aiijj)*sin($(sqΔ/2)*t))+$(anaSol.solPart[1])"
            ss*="\n x2(t)=$(anaSol.coefs[1])*exp($(-α/2)*t)*cos($(sqΔ/2)*t)+$(-anaSol.coefs[2])*exp($(-α/2)*t)*sin($(sqΔ/2)*t)+$(anaSol.solPart[2])"
            ss*="\n return (odeprob,x1,x2) \n end "
            open(path, "a") do io    
                println(io,ss) 
            end
        end
    end
    vectprStr*="]"
    open(path, "a") do io    
        println(io,vectprStr) 
    end
end
############################################################# user space  #####################################################################










 function createProblemBB()
  #=   coefs=[-0.45,-0.56,-0.98,-1.09,-1.28,-1.86,-2.12,-2.78,-3.37,-4.27,-6.32,-8.5,-10.5,-11.6,-13.48,-15.6,-17.05,-19.87,-22.66,0.46,0.55,0.97,1.05,1.22,1.83,2.19,2.75,3.34,4.22,6.32,7.64,8.6,10.6,12.6,14.05,16.87,18.66,21.66,23.66]
   posCoefs=[0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.46,0.55,0.97,1.05,1.22,1.83,2.19,2.75,3.34,4.22,6.32,7.64,8.6,10.6,12.6,14.05,16.87,18.66,21.66,23.66]
   negCoefs=[-0.11,-0.12,-0.13,-0.14,-0.15,-0.16,-0.17,-0.18,-0.19,-0.2,-0.21,-0.22,-0.23,-0.24,-0.25,-0.26,-0.27,-0.28,-0.29,-0.3,-0.31,-0.32,-0.33,-0.34,-0.35,-0.36,-0.37,-0.38,-0.39,-0.45,-0.56,-0.98,-1.09,-1.28,-1.86,-2.12,-2.78,-3.37,-4.27,-6.32,-8.5,-10.5,-11.6,-13.48,-15.6,-17.05,-19.87,-22.66]
   =#
   coefs=[-0.75,-2.0,-12.37,-21.0,0.75,2.0,12.64,27.19]
   posCoefs=[0.2,2.36,7.21,15.4,29.22]
   negCoefs=[-0.956,-2.65,-5.3,-8.09,-12.65,-18.958]
#=    coefs=[-0.75,-12.37,2.0]
   posCoefs=[0.2]
   negCoefs=[-12.65]  =#
   
   #= coefs=[1.0,20.0]
   negCoefs=[-21.0,-1.1] =#
   initConds=[-1.0,-2.0]
   us=1;ul=20
   counter=0
   writeProblemstoFile("AA",coefs,posCoefs,negCoefs,consructJacobiansSystemAA,initConds,us,ul,counter)
end

#createProblemBB()

function createProblemDD()
    #=   coefs=[-0.45,-0.56,-0.98,-1.09,-1.28,-1.86,-2.12,-2.78,-3.37,-4.27,-6.32,-8.5,-10.5,-11.6,-13.48,-15.6,-17.05,-19.87,-22.66,0.46,0.55,0.97,1.05,1.22,1.83,2.19,2.75,3.34,4.22,6.32,7.64,8.6,10.6,12.6,14.05,16.87,18.66,21.66,23.66]
     posCoefs=[0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.46,0.55,0.97,1.05,1.22,1.83,2.19,2.75,3.34,4.22,6.32,7.64,8.6,10.6,12.6,14.05,16.87,18.66,21.66,23.66]
     negCoefs=[-0.11,-0.12,-0.13,-0.14,-0.15,-0.16,-0.17,-0.18,-0.19,-0.2,-0.21,-0.22,-0.23,-0.24,-0.25,-0.26,-0.27,-0.28,-0.29,-0.3,-0.31,-0.32,-0.33,-0.34,-0.35,-0.36,-0.37,-0.38,-0.39,-0.45,-0.56,-0.98,-1.09,-1.28,-1.86,-2.12,-2.78,-3.37,-4.27,-6.32,-8.5,-10.5,-11.6,-13.48,-15.6,-17.05,-19.87,-22.66]
     =#
     coefs=[-0.75,-2.0,-12.37,-21.0,0.75,2.0,12.64,27.19]
     posCoefs=[0.2,2.36,7.21,15.4,29.22]
     negCoefs=[-0.956,-2.65,-5.3,-8.09,-12.65,-18.958]
  #=    coefs=[-0.75,-12.37,2.0]
     posCoefs=[0.2]
     negCoefs=[-12.65]  =#
     
     #= coefs=[1.0,20.0]
     negCoefs=[-21.0,-1.1] =#
     initConds=[-1.0,-2.0]
     us=1;ul=20
     counter=0
     writeComplexProblemstoFile("DD",coefs,posCoefs,negCoefs,consructJacobiansSystemDD,initConds,us,ul,counter)
  end
  
  createProblemDD()