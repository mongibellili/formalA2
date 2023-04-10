
#helper struct that holds dependency data of an event
#= using TimerOutputs
reset_timer!() =#
struct EventDependencyStruct{T,D} #<: AbstractEventDependecy
    id::Int
    evCont::SVector{T,Float64}
    evDisc::SVector{D,Float64}
    evContRHS::SVector{T,Float64}
end

#Non_linear_ordinary_Diff_Equation: to be instanciated by the macro and passed to QSS_data
struct NLODEProblem{T,Z,Y}
    cacheSize::Int
    initConditions::SVector{T,Float64}    
    discreteVars::Vector{Float64}   #discreteVars::MVector{D,Float64}
    #jacobian::SVector{T,SVector{T,Basic}}# jacobian::SVector{T,SVector{T,Float64}} 
    #jacInts::SVector{T,SVector{T,Int}}
    jacInts::Vector{Vector{Int}}
    eqs::Expr#Vector{Expr}
   #=  zceqs::Vector{Expr}
    eventEqus::Vector{Expr} =#
 #=    discreteJacobian::SVector{T,SVector{D,Int}}
    ZC_jacobian::SVector{Z,SVector{T,Int}}
    ZC_jacDiscrete::SVector{Z,SVector{D,Int}} =#
   # discreteJacobian::SVector{T,SVector{D,Basic}}
    #ZC_jacobian::SVector{Z,SVector{T,Basic}}
    zc_SimpleJac::SVector{Z,SVector{T,Int}}

    #ZC_jacDiscrete::SVector{Z,SVector{D,Basic}}
    eventDependencies::SVector{Y,EventDependencyStruct}# 
    #initJac::MVector{T,MVector{T,Float64}}
    initJac :: Vector{Vector{Float64}}
    #SD::SVector{T,SVector{T,Int}}
    SD::Vector{Vector{Int}}
    #= HZ::SVector{Y,SVector{Z,Int}}
    HD::SVector{Y,SVector{T,Int}}
    SZ::SVector{T,SVector{Z,Int}} =#
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
    SD=Vector{Vector{Int}}()
    initJac=[]
    jac = Vector{Vector{SymEngine.Basic}}()
    JacIntVect=Vector{Vector{Int}}()
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
        if @capture(x, y_ = z_)     
           # println("loop")        # x is expre of type lhs=rhs could ve used head == (:=)
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
                       # @timeit "push initial jac" 
                       initJac = Vector{Vector{Float64}}(undef, T)# 
                        for k=1:T
                            push!(jac, zeros(T))# init jac with zeros
                            push!(JacIntVect,[])
                            push!(equs,code0) #init equs with empty expressions
                            push!(SD,[])
                        end
                        
                    else 
                        error("initial conditions already defined! for discrete variables, use the identifier 'discrete' for discrete variables")
                    end
                   # println("init conds")
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
                    extractJac_from_equs(SD,varNum,z,D,usymbols,dsymbols,jac,JacIntVect,jacDiscrete,initJac,discrVars,contVars)
                    ########push!(equs,:($((twoInOneSimplecase(:($(z))))))) # change it to taylor, default of cache given
                    equs[varNum]=:($((twoInOneSimplecase(:($(z))))))
                    #push!(num_cache_equs,1)# to be deleted later
                else #rhs head==call...to be tested later for  math functions and other possible scenarios or user erros
                    #@timeit "changevarname" 
                 #   println("before changenames")
                    z=changeVarNames_to_q_d(z,stateVarName)
                    if length(param)!=0 println("test param");z=replace_parameters(z,param) end
                   # @timeit "extractjac" 
                  # println("before extractjac")
                    extractJac_from_equs(SD,varNum,z,D,usymbols,dsymbols,jac,JacIntVect,jacDiscrete,initJac,discrVars,contVars)                  
                   # push!(num_cache_equs,:($((twoInOne(:($(z),$(cacheSize)))).args[2])))   #number of caches distibuted                                         
                   #@timeit "get0taylor" 
                 #  println("before twoinone")
                   temp=:($((twoInOne(:($(z),1))).args[2]))   #number of caches distibuted   ...no need interpolation and wrap in expr....before was cuz quote....
                    if num_cache_equs<temp 
                          num_cache_equs=temp
                    end 
                   # push!(equs,z)  #twoInone above did change z inside
                    equs[varNum]=z
                    #println("after twoinone")
                end           
            else#end of equations and user enter something weird...handle later
               # error("expression $x: top level contains only expressions 'A=B' or 'if a b' ")#wait until exclude events
            end#end cases inside @capture
       #after capture A=B (init disc equs) we check for 'if statments'
        elseif x isa Expr && x.head==:if   #@capture if did not work
            #syntax: if args[1]  then args[2] else args[3]
            (length(x.args)!=3 && length(x.args)!=2) && error("use format if A>0 B else C or if A>0 B")
            !(x.args[1] isa Expr && x.args[1].head==:call && x.args[1].args[1]==:> && (x.args[1].args[3]==0||x.args[1].args[3]==0.0)) && error("use the format 'if a>0: change if a>b to if a-b>0")
              x.args[1].args[2] isa Number && error("LHS of >  must be be a variable or an expression!")
              x.args[1].args[2]=changeVarNames_to_q_d(x.args[1].args[2],stateVarName)
              extractZCJac_from_equs(x.args[1].args[2],T,D,usymbols,dsymbols,ZCjac,ZCjacDiscrete)
            if x.args[1].args[2].head==:ref  #if one_Var
                #= ifexpr=quote
                            function  $(Symbol(:g_, 2))(q::Vector{Taylor0{Float64}},d::Vector{Float64}, t::Taylor0{Float64},cache::Vector{Taylor0{Float64}})
                                $((twoInOneSimplecase(:($(x.args[1].args[2])))))
                               # return nothing
                              end 
                         end  =#   

                  push!(zcequs,(twoInOneSimplecase(:($(x.args[1].args[2])))))    
                ########################push!(zcequs,:($((twoInOneSimplecase(:($(x.args[1].args[2]))))))) 
              #  push!(num_cache_zcequs,1) #to be deleted later 
            else # if whole expre ops with many vars            
                temp=:($((twoInOne(:($(x.args[1].args[2]),1))).args[2]))   #number of caches distibuted, given 1 place holder for ex.args[2] to be filled inside and returned
                if num_cache_equs<temp 
                      num_cache_equs=temp
                end 
                #= ifexpr=quote
                    function  $(Symbol(:g_, 2))(q::Vector{Taylor0{Float64}},d::Vector{Float64}, t::Taylor0{Float64},cache::Vector{Taylor0{Float64}})
                        $(x.args[1].args[2]) # no need to do twoInone again rhs is already changed inside
                        #return nothing
                    end 
                end     =#
                push!(zcequs,(x.args[1].args[2]))      
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

           #=  temp=:($((twoInOne(:($(x.args[1].args[2]),1))).args[2]))   #number of caches distibuted, given 1 place holder for ex.args[2] to be filled inside and returned
                if num_cache_equs<temp 
                      num_cache_equs=temp
                end  =#
                newPosEventExprToFunc=changeExprToFirstValue(x.args[2].args[1])  #change u[1] to u[1][0]  # pos ev can't be a symbol ...later maybe add check anyway
                push!(eventequs,newPosEventExprToFunc) 
                if x.args[3].args[1] isa Expr 
                    newNegEventExprToFunc=changeExprToFirstValue(x.args[3].args[1])
                    push!(eventequs,newNegEventExprToFunc) 
                else
                    push!(eventequs,x.args[3]) 
                end
            #= eventexpr=quote
                function  $(Symbol(:g_, 2))(q::Vector{Taylor0{Float64}},d::Vector{Float64}, t::Taylor0{Float64},cache::Vector{Taylor0{Float64}})                     
                      #$(x.args[2])
                      $(x.args[2])
                    end
            end =#
            
           #=  eventexpr2=quote
                function  $(Symbol(:g_, 3))(q::Vector{Taylor0{Float64}},d::Vector{Float64}, t::Taylor0{Float64},cache::Vector{Taylor0{Float64}}) 
                     $(x.args[3])                                            
                end
            end =#
              
            
            #after constructing the equations we move to dependencies: we need to change A[n] to An so that they become symbols
            posEvExp = postwalk(a -> a isa Expr && a.head == :ref ? Symbol((a.args[1]), (a.args[2])) : a, x.args[2])
            negEvExp = postwalk(a -> a isa Expr && a.head == :ref ? Symbol((a.args[1]), (a.args[2])) : a, x.args[3])
           #dump(negEvExp)
            indexPosEv = 2 * length(zcequs) - 1 # store events in order
            indexNegEv = 2 * length(zcequs)   
              #------------------pos Event--------------------#
            posEv_disArr= @SVector fill(NaN, D)   #...better than @SVector zeros(D), I can use NaN
            posEv_conArr= @SVector fill(NaN, T)  
            posEv_conArrRHS=@SVector zeros(T)          
            for j = 1:length(posEvExp.args)  # j coressponds the number of statements under one posEvent
                length(posEvExp.args[j].args)!=2 && error("event should be A=B")
                posEvLHS = posEvExp.args[j].args[1]
               !(posEvLHS isa Symbol) && error("lhs of events must be a continuous or a discrete variable")
                
                #!(posEvExp.args[j].args[2] isa Number) #&& error("rhs of events must be a number")#remove the error throw
          


                posEvRHS = posEvExp.args[j].args[2]
                basicRHS = convert(Basic, posEvRHS)
                for l = 1:T             
                    if posEv_conArrRHS[l]==0.0
                        coef = diff(basicRHS, usymbols[l])
                        #posEv_conArrRHS[l]=coef
                        posEv_conArrRHS= setindex(posEv_conArrRHS,coef,l)
                    end
                end




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
            negEv_conArrRHS=@SVector zeros(T)
            if negEvExp.args[1] != :nothing
                #@show negEvExp
                for j = 1:length(negEvExp.args)  # j coressponds the number of statements under one negEvent
                    length(negEvExp.args[j].args)!=2 && error("event should be A=B")
                    negEvLHS = negEvExp.args[j].args[1]
                    !(negEvLHS isa Symbol) && error("lhs of events must be a continuous or a discrete variable")
                    
                    #!(negEvExp.args[j].args[2] isa Number) && error("rhs of events must be a number")
                    negEvRHS = negEvExp.args[j].args[2]
                    basicRHS = convert(Basic, negEvRHS)
                    for l = 1:T             
                        if negEv_conArrRHS[l]==0.0
                            coef = diff(basicRHS, usymbols[l])
                           # negEv_conArrRHS[l]=coef
                            negEv_conArrRHS= setindex(negEv_conArrRHS,coef,l)
                        end
                    end
                    # take care of LHS
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
            structposEvent = EventDependencyStruct(indexPosEv, posEv_conArr, posEv_disArr,posEv_conArrRHS) # right now posEv_conArr is vect of floats
            push!(evsArr, structposEvent)
            structnegEvent = EventDependencyStruct(indexNegEv, negEv_conArr, negEv_disArr,negEv_conArrRHS)
            push!(evsArr, structnegEvent)

        end #end cases inside postwalk
      return x  #
    end #end parent postwalk #########################################################################################################
   # println("end of postwalk")
    Z=length(zcequs)
    Y=2*Z
    if size(jac,1)==T && length(jac[1])==T
        # staticJac = SVector{T,SVector{T,Basic}}(tuple(jac...))
    else
        error("dimension mismatch jac= ",jac," please report the bug")
    end
    if size(jacDiscrete,1)==T && length(jacDiscrete[1])==D
        #staticjacDiscrete=SVector{T,SVector{D,Basic}}(tuple(jacDiscrete...))
    else
        error("dimension mismatch jacDiscrete= ",jacDiscrete," please report the bug")
    end
    
    #println("end of useless convesrion to svector")
    allEpxpr=Expr(:block)
    ##############diffEqua######################
    io = IOBuffer() # i guess this is a sys of diff equ solver so I am not optimizing for case 1 equ
    write(io, "if j==1  $(equs[1]) ;return nothing")
    for i=2:length(equs)
        write(io, " elseif j==$i $(equs[i]) ;return nothing")
    end
    write(io, " end ")
    s = String(take!(io))
    close(io)
    myex1=Meta.parse(s)
    push!(allEpxpr.args,myex1)


     ##############ZCF######################
     
    if length(zcequs)>0
        io = IOBuffer()
        write(io, "if zc==1  $(zcequs[1]) ;return nothing")
        for i=2:length(zcequs)
            write(io, " elseif zc==$i $(zcequs[i]) ;return nothing")
        end
        # write(io, " else  $(zcequs[length(zcequs)]) ;return nothing end ")
        write(io, " end ")
        s = String(take!(io))
        close(io)
        myex2=Meta.parse(s)
        push!(allEpxpr.args,myex2)
    end
    ##############events######################
    if length(eventequs)>0
        io = IOBuffer()
        write(io, "if ev==1  $(eventequs[1]) ;return nothing")
        for i=2:length(eventequs)
            write(io, " elseif ev==$i $(eventequs[i]) ;return nothing")
        end
        # write(io, " else  $(zcequs[length(zcequs)]) ;return nothing end ")
        write(io, " end ")
        s = String(take!(io))
        close(io)
        myex3=Meta.parse(s)
        push!(allEpxpr.args,myex3)
    end

    Base.remove_linenums!(allEpxpr)
    def=Dict{Symbol,Any}()
    def[:head] = :function
    def[:name] = :f   
    def[:args] = [:(j::Int),:(zc::Int),:(ev::Int),:(q::Vector{Taylor0{Float64}}),:(d::Vector{Float64}), :(t::Taylor0{Float64}),:(cache::Vector{Taylor0{Float64}})]
    def[:body] = allEpxpr
    #def[:rtype]=:nothing# test if cache1 always holds sol  
    functioncode=combinedef(def)


   # println("end of function")
  
   # @show jacobian
   # @show d
   #  @show d[1]
  #=   dsym = [symbols("d$i") for i in 1:D]
    qsym = [symbols("q$i") for i in 1:T] =#
   
   #=  qsym = [symbols("q$i") for i in 1:T]
   # @timeit "getinitjac" 
    for i = 1:T
        temparr=Array{Float64}(undef, T)
        for j=1:T  
           ex=jac[i][j]
           for k in eachindex(dsymbols)
               ex=subs(ex, dsymbols[k]=>discrVars[k])#getback d[0] d[1]...in order to get initial correct jacobian to get initial Aii#later check and add sysmbols qi....
           end   
           for k in eachindex(qsym)
            ex=subs(ex, qsym[k]=>contVars[k])#getback d[0] d[1]...in order to get initial correct jacobian to get initial Aii#later check and add sysmbols qi....
           end    
            temparr[j]=ex
        end
        tempJac[i]=temparr
        
     end =#
    # println("end of substitution")
    # initJac = MVector{T,MVector{T,Float64}}(tuple(tempJac...))

     #println("create dependencies")
    # do i really need to use staticarrays here, since i am going to  use them once (below)?????
    staticZC_jacobian=SVector{Z,SVector{T,Basic}}(tuple(ZCjac...))
    staticZC_jacDiscrete=SVector{Z,SVector{D,Basic}}(tuple(ZCjacDiscrete...))
    eventDependencies = SVector{Y,EventDependencyStruct}(evsArr)  # 2*Z each zc yields 2 events
   
      #usedto track if zc changed sign; each zc has a value and a sign 
     # @timeit "changejacbasicToInts" jacInts=changeBasicToInts(staticJac)# change from type nonisbits to int so that access is cheaper down
     zc_SimpleJac=changeBasicToInts(staticZC_jacobian)
     #*******************************create dependencies**************************
    # @timeit "createSD" SD = createDependencyMatrix(staticJac)
    #=  dD =createDependencyMatrix(jacDiscrete) # temp dependency to be used to determine HD1 and HZ1 HD=Hd-dD Union Hs-sD
    # @timeit "createSZ"
      SZ =createDependencyMatrix(staticZC_jacobian) 
     dZ =createDependencyMatrix(staticZC_jacDiscrete) =# # temp dependency to be used to determine HD2 and HZ2
    # HZ1HD1=createDependencyToEventsDiscr(dD,dZ,eventDependencies) 
    # HZ2HD2=createDependencyToEventsCont(SD,SZ,eventDependencies) 
    # HZ=unionDependency(HZ1HD1[1],HZ2HD2[1])
    # @timeit "uniondependecy" HD=unionDependency(HZ1HD1[2],HZ2HD2[2])initJac
    #@show functioncode
   # print_timer()
  # println("build struct")
    myodeProblem = NLODEProblem(num_cache_equs,contVars, discrVars, JacIntVect,functioncode, zc_SimpleJac, eventDependencies,initJac,SD#= ,HZ,HD,SZ =#)
   # path="/home/unknown/relaxedqssA/relaxedqssA/src/models/loopToBody.jl"
    
   # open(path, "a") do io  println(io,string(myodeProblem))  end
  # (functioncode,myodeProblem )
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