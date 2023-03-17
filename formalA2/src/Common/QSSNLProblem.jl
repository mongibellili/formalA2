
#helper struct that holds dependency data of an event
struct EventDependencyStruct{T,D} #<: AbstractEventDependecy
    id::Int
    evCont::SVector{T,Float64}
    evDisc::SVector{D,Float64}
end

#Non_linear_ordinary_Diff_Equation: to be instanciated by the macro and passed to QSS_data
struct NLODEProblem{T,D,Z,Y}
    cacheSize::Int
    initConditions::SVector{T,Float64}    
    discreteVars::Vector{Float64}   #discreteVars::MVector{D,Float64}
    jacobian::SVector{T,SVector{T,Basic}}# jacobian::SVector{T,SVector{T,Float64}} 
    eqs::Expr#Vector{Expr}
    zceqs::Vector{Expr}
    eventEqus::Vector{Expr}
    discreteJacobian::SVector{T,SVector{D,Basic}}
    ZC_jacobian::SVector{Z,SVector{T,Basic}}
    ZC_jacDiscrete::SVector{Z,SVector{D,Basic}}
    eventDependencies::SVector{Y,EventDependencyStruct}# 
end
#macro receives user code and creates the problem struct
function NLodeProblem(odeExprs)
    #Base.remove_linenums!(odeExprs)
    stateVarName=:u #default
    du=:du #default
    initCondreceived=false #bool throw error if user redefine
    #discreteVarName=:d
    T=0#numberStateVars=
    D=0#numberDiscreteVars=
    Z=0#numberzcFunctions=
    contVars=[]
    discrVars=[]
    jac = Vector{Vector{SymEngine.Basic}}()
    jacDiscrete = Vector{Vector{SymEngine.Basic}}()
    ZCjac = Vector{Vector{SymEngine.Basic}}()
    ZCjacDiscrete = Vector{Vector{SymEngine.Basic}}()
    usymbols=[]
    xsymbols=[]#added in case later zcf(x,d,t) is better than zcf(q,d,t)
    dsymbols=[]
    param=Dict{Symbol,Number}()
    equs=Vector{Expr}()# vector to collect diff-equations
    num_cache_equs=1#cachesize
    zcequs=Vector{Expr}()#vect to collect if-statements
    eventequs=Vector{Expr}()#vect to collect events
    evsArr = [] #helper to collect info about event dependencies
    postwalk(odeExprs) do x
        if @capture(x, y_ = z_)             # x is expre of type lhs=rhs could ve used head == (:=)
            if y isa Symbol && z isa Number                    
                param[y]=z # parameters              
            elseif z isa Expr && z.head==:vect # rhs ==vector of state vars initCond or discrete vars
                if y!=:discrete
                    if !initCondreceived
                        initCondreceived=true
                        stateVarName=y
                        du=Symbol(:d,stateVarName)
                        T = length(z.args)
                        contVars = SVector{T,Float64}(z.args) 
                        usymbols = [symbols("q$i") for i in 1:T] # symbols for cont vars
                        xsymbols = [symbols("x$i") for i in 1:T] # symbols for cont vars
                        code0=Expr(:block)
                        for k=1:T
                            push!(jac, zeros(T))# init jac with zeros
                            push!(equs,code0) #init equs with empty expressions
                        end
                        
                    else 
                        error("initial conditions already defined! for discrete variables, use the identifier 'discrete' for discrete variables")
                    end
                else #y==:discrete #later forbid redefine discrete like cont
                    D = length(z.args)
                    discrVars = Vector{Float64}(z.args)#discrVars = MVector{D,Float64}(z.args)
                    dsymbols = [symbols("d$i") for i in 1:D] # symbols for disc vars
                    for k=1:T
                        push!(jacDiscrete, zeros(D)) 
                    end

                end
            elseif  du==y.args[1] && ( (z isa Expr && (z.head==:call || z.head==:ref)) || z isa Number)#z is rhs of diff-equations because du==
                varNum=y.args[2] # order of variable
               
               
               
               if z isa Number # rhs of equ =number  
                    #= push!(jac, zeros(T)) #no dependencies
                    push!(jacDiscrete, zeros(D))  =#
                   ###### push!(equs,:($((twoInOneSimplecase(:($(z))))))) # change it to taylor
                   equs[varNum]=:($((twoInOneSimplecase(:($(z))))))
                   # push!(num_cache_equs,1) #was thinking about internal cache_clean_up...hurt performance...to be deleted later
                elseif z.head==:ref #rhs is only one var
                    z=changeVarNames_to_q_d(z,stateVarName)# the user may choose any symbols for continuos only, discrete naming is fixed to eliminate ambiguities
                    extractJac_from_equs(varNum,z,T,D,usymbols,dsymbols,jac,jacDiscrete)
                    ########push!(equs,:($((twoInOneSimplecase(:($(z))))))) # change it to taylor, default of cache given
                    equs[varNum]=:($((twoInOneSimplecase(:($(z))))))
                    #push!(num_cache_equs,1)# to be deleted later
                else #rhs head==call...to be tested later for  math functions and other possible scenarios or user erros
                    z=changeVarNames_to_q_d(z,stateVarName)
                    z=replace_parameters(z,param)
                    extractJac_from_equs(varNum,z,T,D,usymbols,dsymbols,jac,jacDiscrete)                  
                   # push!(num_cache_equs,:($((twoInOne(:($(z),$(cacheSize)))).args[2])))   #number of caches distibuted                                         
                    temp=:($((twoInOne(:($(z),1))).args[2]))   #number of caches distibuted   ...no need interpolation and wrap in expr....before was cuz quote....
                    if num_cache_equs<temp 
                          num_cache_equs=temp
                    end 
                   # push!(equs,z)  #twoInone above did change z inside
                    equs[varNum]=z
                end           
            else#end of equations and user enter something weird...handle later
               # error("expression $x: top level contains only expressions 'A=B' or 'if a b' ")#wait until exclude events
            end#end cases inside @capture
       #after capture A=B (init disc equs) we check for 'if statments'
        elseif x isa Expr && x.head==:if   #@capture if did not work
            #syntax: if args[1]  then args[2] else args[3]
            (length(x.args)!=3 && length(x.args)!=2) && error("use format if A>0 B else C")
            !(x.args[1] isa Expr && x.args[1].head==:call && x.args[1].args[1]==:> && (x.args[1].args[3]==0||x.args[1].args[3]==0.0)) && error("use the format 'if a>0")
              x.args[1].args[2] isa Number && error("LHS of term must be be a variable or an expression!")
              x.args[1].args[2]=changeVarNames_to_q_d(x.args[1].args[2],stateVarName)
              extractZCJac_from_equs(x.args[1].args[2],T,D,usymbols,dsymbols,ZCjac,ZCjacDiscrete)
            if x.args[1].args[2].head==:ref  #if one_Var
                ifexpr=quote
                            function  $(Symbol(:g_, 2))(q::Vector{Taylor0{Float64}},d::Vector{Float64}, t::Taylor0{Float64},cache::Vector{Taylor0{Float64}})
                                $((twoInOneSimplecase(:($(x.args[1].args[2])))))
                              end 
                         end    
                  push!(zcequs,ifexpr)     
                ########################push!(zcequs,:($((twoInOneSimplecase(:($(x.args[1].args[2]))))))) 
              #  push!(num_cache_zcequs,1) #to be deleted later 
            else # if whole expre ops with many vars            
                temp=:($((twoInOne(:($(x.args[1].args[2]),1))).args[2]))   #number of caches distibuted, given 1 place holder for ex.args[2] to be filled inside and returned
                if num_cache_equs<temp 
                      num_cache_equs=temp
                end 
                ifexpr=quote
                    function  $(Symbol(:g_, 2))(q::Vector{Taylor0{Float64}},d::Vector{Float64}, t::Taylor0{Float64},cache::Vector{Taylor0{Float64}})
                        $(x.args[1].args[2]) # no need to do twoInone again rhs is already changed inside
                    end 
                end    
                push!(zcequs,ifexpr)      
            end     
            ##################################################################################################################
            #                                                      events     
            ##################################################################################################################                  
            # each 'if-statmets' has 2 events (arg[2]=posEv and arg[3]=NegEv) each pos or neg event has a function...later i can try one event for zc
            x.args[2]=changeVarNames_to_q_d(x.args[2],stateVarName) #posEvent

            if length(x.args)==2  #if user only wrote the positive evnt, here I added the negative event wich does nothing
               
                equalexpr = quote nothing end # neg dummy event 
                push!(x.args, equalexpr)
                Base.remove_linenums!(x.args[3])
            else
               x.args[3]=changeVarNames_to_q_d(x.args[3],stateVarName) #negEvent
            end
            eventexpr=quote
                function  $(Symbol(:g_, 2))(q::Vector{Taylor0{Float64}},d::Vector{Float64}, t::Taylor0{Float64},cache::Vector{Taylor0{Float64}})                     
                      $(x.args[2])
                    end
            end
            push!(eventequs,eventexpr)  
            eventexpr2=quote
                function  $(Symbol(:g_, 3))(q::Vector{Taylor0{Float64}},d::Vector{Float64}, t::Taylor0{Float64},cache::Vector{Taylor0{Float64}}) 
                     $(x.args[3])                                            
                end
            end
            push!(eventequs,eventexpr2)   
              
            #after constructing the equations we move to dependencies: we need to change A[n] to An so that they become symbols
            posEvExp = postwalk(a -> a isa Expr && a.head == :ref ? Symbol((a.args[1]), (a.args[2])) : a, x.args[2])
            negEvExp = postwalk(a -> a isa Expr && a.head == :ref ? Symbol((a.args[1]), (a.args[2])) : a, x.args[3])
#dump(negEvExp)
            indexPosEv = 2 * length(zcequs) - 1 # store events in order
            indexNegEv = 2 * length(zcequs)   
              #------------------pos Event--------------------#
            posEv_disArr= @SVector fill(NaN, D)   #...better than @SVector zeros(D), I can use NaN
            posEv_conArr= @SVector fill(NaN, T)            
            for j = 1:length(posEvExp.args)  # j coressponds the number of statements under one posEvent
                length(posEvExp.args[j].args)!=2 && error("event should be A=B")
               !(posEvExp.args[j].args[1] isa Symbol) && error("lhs of events must be a continuous or a discrete variable")
                posEvLHS = posEvExp.args[j].args[1]
                #!(posEvExp.args[j].args[2] isa Number) #&& error("rhs of events must be a number")#remove the error throw
                posEvRHS = posEvExp.args[j].args[2] #
                basicLHS = convert(Basic, posEvLHS) # because usymbols and dsymbols are of type Basic...also used elsewhere to find derivatives...maybe later throw eroor?
                discVarpositionArray = indexin(basicLHS, dsymbols)#basicLHS is a symbol, dsymbols is a vect of symbols=[d1,d2,d3]    #later try findall(x->x == basicLHS, d)
                #indexin(a, b) Returns a vector containing the highest index in b for each value in a that is a member of b         
                if !(discVarpositionArray[1] === nothing)
                    posEv_disArr = setindex(posEv_disArr, 1.0, discVarpositionArray[1])
                else # lhs is not a disc var 
                    conVarpositionArray = indexin(basicLHS, usymbols)
                    if !(conVarpositionArray[1] === nothing)
                        posEv_conArr = setindex(posEv_conArr, 1.0, conVarpositionArray[1])
                    else
                        error("LHS is neither a cont nor a discr var!!") #later throw error of used variable not declared in initcond vector or discrete vector
                    end
                end
            end

            #------------------neg Event--------------------#
            negEv_disArr= @SVector fill(NaN, D)
            negEv_conArr=@SVector fill(NaN, T)
            if negEvExp.args[1] != :nothing
                #@show negEvExp
                for j = 1:length(negEvExp.args)  # j coressponds the number of statements under one negEvent
                    length(negEvExp.args[j].args)!=2 && error("event should be A=B")
                    !(negEvExp.args[j].args[1] isa Symbol) && error("lhs of events must be a continuous or a discrete variable")
                    negEvLHS = negEvExp.args[j].args[1]
                    #!(negEvExp.args[j].args[2] isa Number) && error("rhs of events must be a number")
                    negEvRHS = negEvExp.args[j].args[2]
                    basicLHS = convert(Basic, negEvLHS)
                    discVarpositionArray = indexin(basicLHS, dsymbols)
                    if !(discVarpositionArray[1] === nothing)
                        negEv_disArr = setindex(negEv_disArr, 1.0, discVarpositionArray[1])
                    else # lhs is not a disc var 
                        conVarpositionArray = indexin(basicLHS, usymbols)
                        if !(conVarpositionArray[1] === nothing)
                            negEv_conArr = setindex(negEv_conArr, 1.0, conVarpositionArray[1])
                        else
                            println("LHS is neither a continous nor a discrete variable!!")
                        end
                    end
                end
            

            end
            structposEvent = EventDependencyStruct(indexPosEv, posEv_conArr, posEv_disArr) # right now posEv_conArr is vect of floats
            push!(evsArr, structposEvent)
            structnegEvent = EventDependencyStruct(indexNegEv, negEv_conArr, negEv_disArr)
            push!(evsArr, structnegEvent)

        end #end cases inside postwalk
      return x  #
    end #end parent postwalk 

    Z=length(zcequs)
    Y=2*Z
    if size(jac,1)==T && length(jac[1])==T
         staticJac = SVector{T,SVector{T,Basic}}(tuple(jac...))
    else
        error("dimension mismatch jac= ",jac," please report the bug")
    end
    if size(jacDiscrete,1)==T && length(jacDiscrete[1])==D
        staticjacDiscrete=SVector{T,SVector{D,Basic}}(tuple(jacDiscrete...))
    else
        error("dimension mismatch jacDiscrete= ",jacDiscrete," please report the bug")
    end
    staticZC_jacobian=SVector{Z,SVector{T,Basic}}(tuple(ZCjac...))
    staticZC_jacDiscrete=SVector{Z,SVector{D,Basic}}(tuple(ZCjacDiscrete...))
    eventDependencies = SVector{Y,EventDependencyStruct}(evsArr)  # 2*Z each zc yields 2 events
    io = IOBuffer() # i guess this is a sys of diff equ solver so I am not optimizing for case 1 equ
    write(io, "if j==1  $(equs[1]) ;return nothing")
    for i=2:length(equs)-1
        write(io, " elseif j==$i $(equs[i]) ;return nothing")
    end
    write(io, " else  $(equs[length(equs)]) ;return nothing end ")
    s = String(take!(io))
    close(io)
    def=Dict{Symbol,Any}()
    def[:head] = :function
    def[:name] = :f   
    def[:args] = [:(j::Int),:(q::Vector{Taylor0{Float64}}),:(d::Vector{Float64}), :(t::Taylor0{Float64}),:(cache::Vector{Taylor0{Float64}})]
    def[:body] = Meta.parse(s)
    #def[:rtype]=:nothing# test if cache1 always holds sol  
    functioncode=combinedef(def)
    myodeProblem = NLODEProblem(num_cache_equs,contVars, discrVars, staticJac, functioncode,zcequs,eventequs, staticjacDiscrete, staticZC_jacobian, staticZC_jacDiscrete, eventDependencies)
   # path="/home/unknown/relaxedqssA/relaxedqssA/src/models/loopToBody.jl"
    
   # open(path, "a") do io  println(io,string(myodeProblem))  end
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
   
    #path="/home/unknown/relaxedqssA/relaxedqssA/src/models/loopToBody.jl"
    path="./ModelsDiffEq.jl"
    def[:body] = code;newFun=combinedef(def)
    open(path, "a") do io  println(io,string(newFun))  end
  end