using MacroTools: postwalk, prewalk, @capture
using ExprTools
function NLodeProblem(odeExprs)
   dump(odeExprs)
   @show odeExprs

end

macro NLodeProblem(odeExprs)
    Base.remove_linenums!(odeExprs)
    equs=Vector{Expr}()
    code0=Expr(:block)
    for i=1:2
        push!(equs,code0)
    end
    equs[1]=odeExprs
    #@show equs
    
    #dump(odeExprs)

#=    indexArgsLoop=0
    postwalk(odeExprs) do loop
        if @capture(loop, for var_ in 1:niter_ loopbody__ end)
            indexArgsLoop = findall(x->x==loop, odeExprs.args)[1]
            #indexArgsLoop= indexin(loop, odeExprs.args)#basicLHS is a symbol, dsymbols is a vect of symbols=[d1,d2,d3]    #later try findall(x->x == basicLHS, d)
            #indexin(a, b) Returns a vector containing the highest index in b for each value in a that is a member of b         
           @show indexArgsLoop
        end
        return loop
    end
    if indexArgsLoop!=0
    sizeargs=length(odeExprs.args)
     temp=odeExprs.args[sizeargs]
    odeExprs.args[sizeargs]=odeExprs.args[indexArgsLoop]
    odeExprs.args[indexArgsLoop]=temp
    pop!(odeExprs.args)
    end =#
    # dump(odeExprs)
    #=     postwalk(odeExprs) do loop
        if @capture(loop, for var_ in 1:niter_ loopbody__ end)
           # @show niter
            for i in 1:niter
                ex = postwalk(a -> a  == var ? i : a, loopbody[1]) # obtain ex = the body of the loop with counter changed to actual number
                v=postwalk(ex) do a  # this is to change for example 5+1 to 6 in vect indices
                    if a isa Expr && a.head == :ref && a.args[2] isa Expr
                        a.args[2]=eval(a.args[2])
                    end
                    return a
                end#end interior postwalk
              #  Base.remove_linenums!(v)# one equation v to be cleaned then saved (pushed alone to the parent expression)
                push!(odeExprs.args,v)
            end#end for loop
        end#end if capture
      return loop
    end#end postwalk look for "for loop"
    # code = Expr(:block)
    #=  code=quote fm=5 end
    # @show odeExprs
    posEvExp = postwalk(a -> a isa Expr && a.head == :for ? code.args[2] : a, odeExprs) =#
    newex=postwalk(odeExprs) do x
        if @capture(x, for var_ in 1:niter_ loopbody__ end)       
                    x=postwalk(x) do y
                        
                        if y isa Expr 
                            
                            y.head=:block
                        end
                        return y
                    end
                   # dump(x)
        end#end if capture  
        return x
    end
    # dump(posEvExp)
    NLodeProblem(newex) =#


   # 
  NLodeProblem(equs)


end

#= @loopToBody function foo()
    for i in 1:2  
         du[i]=(((U/Rs[i]) - iL[i]) * (Rs[i]*Rd[i]/(Rs[i]+Rd[i])) - uC)/L;
    end         
end =#
@NLodeProblem begin
   
        parameter2=0.00001 ;
    
    #=  u = [1.0, -1.3,0.5]
        discrete = [0.0] =#
       #=  du[1] = -u[1]+u[2]
        du[2] =u[3]-u[2]   =# 
        du[3]=-u[1]
        for i in 1:2  
            du[i]=u[i+1]-u[i]
       end 
      
   
     #=    parameter1=3000.0# cache can be dynamic....parameters take this feature
        parameter2=0.00001 =#
 #=  if u[1]>0   #5*discrte gave error
    discrete[1]=0.0   #discrete=0.0-->type Symbol has no field args...find to personalize error msg            
else
    discrete[1]=1.0                                    
end =#

end