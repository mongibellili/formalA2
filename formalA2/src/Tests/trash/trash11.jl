using ExprTools
allEpxpr=Expr(:block)
#allEpxpr=:()
io = IOBuffer() # i guess this is a sys of diff equ solver so I am not optimizing for case 1 equ
    ##############diffEqua######################
    write(io, "if j==1  case1 ;return nothing")
    for i=2:3
        write(io, " elseif j==$i case$i ;return nothing")
    end
   # write(io, " else  $(equs[length(equs)]) ;return nothing end ")
    write(io, " end ")
     

    s = String(take!(io))
    close(io)

    myex1=Meta.parse(s)
   # @show myex1
    push!(allEpxpr.args,myex1)
    io = IOBuffer() # i guess this is a sys of diff equ solver so I am not optimizing for case 1 equ
 
     ##############ZCF######################
     #write(io, "\n")
    write(io, "if zc==1  zscase1 ;return nothing")
    for i=2:3
        write(io, " elseif zc==$i zccase$i ;return nothing")
    end
   # write(io, " else  $(zcequs[length(zcequs)]) ;return nothing end ")
   write(io, " end ")

    s = String(take!(io))
    close(io)

    myex2=Meta.parse(s)
    #@show myex2
    push!(allEpxpr.args,myex2)
    Base.remove_linenums!(allEpxpr)
    #@show allEpxpr

    def=Dict{Symbol,Any}()
    def[:head] = :function
    def[:name] = :f1
    def[:args] = [:(u::Vector{Float64})]
    def[:body] = allEpxpr
    

    path="./ModelsDiffEq33.jl"
    newFun=combinedef(def)
    open(path, "a") do io  println(io,string(newFun))  end