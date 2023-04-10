
using SymEngine
using StaticArrays
macro test()
m=:((1.0 - 2*q1 - q2)*100.8)
ex = convert(Basic, m)
newex=Array{Float64}(undef, 2)
contVars = [1.2,2.0]
xsymbols = [symbols("x$i") for i in 1:2]
for k in eachindex(xsymbols)
    display(ex)
    ex=subs(ex, xsymbols[k]=>contVars[k])#getback d[0] d[1]...in order to get initial correct jacobian to get initial Aii#later check and add sysmbols qi....
end    
newex[1]=ex
display(typeof(newex))
@show newex
end
@test()