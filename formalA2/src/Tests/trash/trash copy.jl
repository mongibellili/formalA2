#using MacroTools
using MacroTools: postwalk, prewalk, @capture
using ExprTools

macro loopToBody(loopfun)
  Base.remove_linenums!(loopfun)
  loopfun.head != :function && error("Expression is not a function definition!")
  def=splitdef(loopfun)
  loop=def[:body]
  @assert(@capture(loop, for var_ in 1:niter_ loopbody__ end),"macro for loop")
  code = Expr(:block)
  #dump(loopbody[1];maxdepth=16)
  
#@show loopbody[1]
  for i in 1:niter
    #v=:(a=5)
    #v=:($loopbody)
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
 
  path="/home/unknown/relaxedqssA/relaxedqssA/src/models/loopToBody.jl";def[:body] = code;newFun=combinedef(def)
  open(path, "a") do io  println(io,string(newFun))  end
end
#a=zeros(5)
@loopToBody function foo()
    for k in 1:5   
         du[k]=3*k*u[k+1]-2*u[k]
    end         
end
#display(prob) 