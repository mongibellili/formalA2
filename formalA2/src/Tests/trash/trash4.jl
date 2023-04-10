using InteractiveUtils
#= u=[[2.0,3.33,5.6],[12.111]]
v=[12.0,12.8,3.23,6.0001] =#
#= for b in (u[1]  )  
    @show b
end =#

function fun1()
    u=[[2.0,3.33,5.6],[12.111]]
    for i in u[1]
        #@show 
        i
    end
end
function fun2()
    u=[[2.0,3.33,5.6],[12.111]]
    for l = 1:length(u[1])
        #@show 
        u[1][l]
    end
end
fun1()
fun2()
display(@code_warntype fun2()) 