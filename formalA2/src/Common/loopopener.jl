#= using MacroTools: postwalk, prewalk, @capture
using ExprTools
function problemODE(code)
   dump(code)
end

macro loopToBody(loopfun)
  Base.remove_linenums!(loopfun)
  loopfun.head != :function && error("Expression is not a function definition!")
  def=splitdef(loopfun)
  loop=def[:body]
  @assert(@capture(loop, for var_ in 1:niter_ loopbody__ end),"macro for loop")
  code = Expr(:block)
  for i in 1:niter
    ex = postwalk(a -> a  == var ? i : a, loopbody[1])
    v=postwalk(ex) do a  # change rhs equ to contain q[] and d[] instead of user symbols
        if a isa Expr && a.head == :ref && a.args[2] isa Expr
            a.args[2]=eval(a.args[2])
        end
        return a
    end
    Base.remove_linenums!(v)
    push!(code.args,v)
  end
 problemODE(code)
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
@loopToBody function foo()
  for i in 1:2  
       du[i]=u[i+1]-u[i];
  end         
end =#