using MacroTools: postwalk, prewalk, @capture
using ExprTools
function problemODE(code)
   dump(code)
end

macro odeproblem(odeExprs)
  Base.remove_linenums!(odeExprs)
  postwalk(odeExprs) do x
    if @capture(x, for var_ in 1:niter_ loopbody__ end)
   # code = Expr(:block)
   @show niter
        #= for i in 1:niter
            ex = postwalk(a -> a  == var ? i : a, loopbody[1]) # obtain ex = the body of the loop with counter changed to actual number
            v=postwalk(ex) do a  # this is to change for example 5+1 to 6 in vect indices
                if a isa Expr && a.head == :ref && a.args[2] isa Expr
                    a.args[2]=eval(a.args[2])
                end
                return a
            end#end interior postwalk
            Base.remove_linenums!(v)# one equation v to be cleaned then saved (pushed alone to the parent expression)
            # push!(code.args,v)
             push!(odeExprs.args,v)
        end#end for loop =#
    end#end if capture
  end#end postwalk look for "for loop"

 #problemODE(odeExprs)
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
@odeproblem begin
    parameter1=3000.0# cache can be dynamic....parameters take this feature
        parameter2=0.00001
        u = [10.0,0.0]
  for i in 1:2  
       du[i]=u[i+1]-u[i];
  end         
end