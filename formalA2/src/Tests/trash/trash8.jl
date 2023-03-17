using MacroTools: postwalk, prewalk, @capture
using ExprTools
function NLodeProblem(odeExprs)
   #dump(odeExprs)
   @show odeExprs

end

macro NLodeProblem(odeExprs)
    Base.remove_linenums!(odeExprs)
    code = Expr(:block)
    postwalk(odeExprs) do loop
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
                push!(code.args,v)
            end#end for loop
        end#end if capture
      return loop
    end#end postwalk look for "for loop"
   # code = Expr(:block)
   for i=1:length(code.args)
    push!(odeExprs.args,code.args[i])
   end
   #=  code=quote fm=5 end
   # @show odeExprs
    posEvExp = postwalk(a -> a isa Expr && a.head == :for ? code.args[2] : a, odeExprs) =#
   #=  indexArgsLoop=0
    postwalk(odeExprs) do loop
        if @capture(loop, for var_ in 1:niter_ loopbody__ end)
            indexArgsLoop = findall(x->x==loop, odeExprs.args)[1]
            #indexArgsLoop= indexin(loop, odeExprs.args)#basicLHS is a symbol, dsymbols is a vect of symbols=[d1,d2,d3]    #later try findall(x->x == basicLHS, d)
            #indexin(a, b) Returns a vector containing the highest index in b for each value in a that is a member of b         
           #@show indexArgsLoop
        end
        return loop
    end =#
    
   #=  if indexArgsLoop!=0 # if for loop exists
        sizeargs=length(odeExprs.args)
        if indexArgsLoop==sizeargs #if for loop is the last arg
            pop!(odeExprs.args)
        else

            temp=odeExprs.args[sizeargs]#backup last arg
            odeExprs.args[sizeargs]=odeExprs.args[indexArgsLoop]# put for loop at end
            odeExprs.args[indexArgsLoop]=temp# item backed up goes to index(place) of for loop
            pop!(odeExprs.args)# remove for loop
        end
    end =#
  # dump(posEvExp)
  NLodeProblem(odeExprs)
  #path="/home/unknown/relaxedqssA/relaxedqssA/src/models/loopToBody.jl"
#=   path="./loopToBody.jl"
  def[:body] = code;newFun=combinedef(def)
  open(path, "a") do io  println(io,string(newFun))  end =#
end

#= @loopToBody function foo()
    for i in 1:2  
         du[i]=(((U/Rs[i]) - iL[i]) * (Rs[i]*Rd[i]/(Rs[i]+Rd[i])) - uC)/L;
    end         
end =#
@NLodeProblem begin
    
       #=  u = [1.0, -1.3,0.5]
        discrete = [0.0] =#
       #=  du[1] = -u[1]+u[2]
        du[2] =u[3]-u[2]   =# 
     
        du[3]=-u[1]
        parameter1=3000.0# cache can be dynamic....parameters take this feature
        parameter2=0.00001
        for i in 1:2  
            du[i]=u[i+1]-u[i];
       end 
 #=  if u[1]>0   #5*discrte gave error
    discrete[1]=0.0   #discrete=0.0-->type Symbol has no field args...find to personalize error msg            
else
    discrete[1]=1.0                                    
end =#

end