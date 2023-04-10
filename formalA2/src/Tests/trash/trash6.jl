
using TimerOutputs
reset_timer!()
 macro twostep(arg)
    u=[]
   @timeit "push" push!(u,2.33)
   @timeit "print" println("I execute at parse time. The argument is: ", u)
    return :(println("I execute at runtime. The argument is: ", $arg))
end


macro test()
@twostep(1)
println("hi")
end

print_timer()

#@test()
#@twostep(2)