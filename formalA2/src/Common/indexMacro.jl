macro NLodeProblem(odeExprs)
    Base.remove_linenums!(odeExprs)
   # println("begin normal macro")
    NLodeProblem(odeExprs)
  #path="/home/unknown/relaxedqssA/relaxedqssA/src/models/loopToBody.jl"
#=   path="./loopToBody.jl"
  def[:body] = code;newFun=combinedef(def)
  open(path, "a") do io  println(io,string(newFun))  end =#
end
macro NLodeProblemLoop(odeExprs)
   # println("begin loop unrol")
    Base.remove_linenums!(odeExprs)
    code = Expr(:block)
    postwalk(odeExprs) do loop
        if @capture(loop, for var_ in b_:niter_ loopbody__ end)
            for i in b:niter
                ex = postwalk(a -> a  == var ? i : a, loopbody[1]) # obtain ex = the body of the loop with counter changed to actual number
                v=postwalk(ex) do a  # this is to change for example 5+1 to 6 in vect indices
                    if a isa Expr && a.head == :ref && a.args[2] isa Expr # counter outside ref still not handeled...init cond or used in equations.later***************************************************************
                        a.args[2]=eval(a.args[2])
                    end
                    return a
                end#end interior postwalk
                push!(code.args,v)# must not change odeexprs args array here : memory crash
            end#end for loop
        end#end if capture
     return loop
    end#end postwalk look for "for loop"
    for i=1:length(code.args)# used another var to store the equation: 
        push!(odeExprs.args,code.args[i])
    end
   # section that deletes not needed for loop
    indexArgsLoop=0
    postwalk(odeExprs) do loop
        if @capture(loop, for var_ in b_:niter_ loopbody__ end)
            indexArgsLoop = findall(x->x==loop, odeExprs.args)[1]
        end
        return loop
    end
    if indexArgsLoop!=0
    sizeargs=length(odeExprs.args)
     temp=odeExprs.args[sizeargs]
    odeExprs.args[sizeargs]=odeExprs.args[indexArgsLoop]
    odeExprs.args[indexArgsLoop]=temp
    pop!(odeExprs.args)
    end
    #e,nd of section that deletes for loop
   # @show odeExprs
  #println("end unroll loop")
  NLodeProblem(odeExprs)
  #path="/home/unknown/relaxedqssA/relaxedqssA/src/models/loopToBody.jl"
#=   path="./loopToBody.jl"
  def[:body] = code;newFun=combinedef(def)
  open(path, "a") do io  println(io,string(newFun))  end =#
end