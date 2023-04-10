#using TimerOutputs
function QSS_integrate(::Val{O}, s::QSS_data{T,Z}, odep::NLODEProblem{T,Z,Y},f::Function) where {O,T,Z,Y}
  
  #reset_timer!()
  ft = s.finalTime;initTime = s.initialTime;relQ = s.dQrel;absQ = s.dQmin;maxErr=s.maxErr;savetimeincrement=s.savetimeincrement;savetime = savetimeincrement
  #*********************************qss method data*****************************************
  quantum = s.quantum;nextStateTime = s.nextStateTime;nextEventTime = s.nextEventTime;nextInputTime = s.nextInputTime
  tx = s.tx;tq = s.tq;x = s.x;q = s.q;t=s.t
  savedVars=s.savedVars;savedTimes=s.savedTimes;integratorCache=s.integratorCache;taylorOpsCache=s.taylorOpsCache;cacheSize=odep.cacheSize
  a=odep.initJac;
  #a=odep.tempJac;
 
 # u=s.u;tu=s.tu
  #*********************************problem info*****************************************
  d = odep.discreteVars
  jac=odep.jacInts
  #display(jac);println()
 #=  eqs::Expr#Vector{Expr}
  zceqs::Vector{Expr}
  eventEqus::Vector{Expr} =#
#=    discreteJacobian::SVector{T,SVector{D,Int}}
  ZC_jacobian::SVector{Z,SVector{T,Int}}
  ZC_jacDiscrete::SVector{Z,SVector{D,Int}} =#
 # discreteJacobian::SVector{T,SVector{D,Basic}}
  #ZC_jacobian::SVector{Z,SVector{T,Basic}}
  zc_SimpleJac=odep.zc_SimpleJac

  #ZC_jacDiscrete::SVector{Z,SVector{D,Basic}}

 
  SD=odep.SD
  @show SD
 #=  HZ=odep.HZ
  HD=odep.HD
  SZ=odep.SZ =#
 # jacobian = odep.jacobian
  #discJac = odep.discreteJacobian
 # zc_jac = odep.ZC_jacobian
  #ZC_discJac = odep.ZC_jacDiscrete
  evDep = odep.eventDependencies
  #********************************helper values*******************************  
 # qaux=s.qaux;olddx=s.olddx;olddxSpec = zeros(MVector{T,MVector{O,Float64}}) # later can only care about 1st der
  numSteps = zeros(MVector{T,Int})
  oldsignValue = MMatrix{Z,2}(zeros(Z*2))  #usedto track if zc changed sign; each zc has a value and a sign 
  #= jac=changeBasicToInts(jacobian)# change from type nonisbits to int so that access is cheaper down
  zc_SimpleJac=changeBasicToInts(zc_jac) =#
  #*******************************create dependencies**************************
#=   SD = createDependencyMatrix(jacobian)
  dD =createDependencyMatrix(discJac) # temp dependency to be used to determine HD1 and HZ1 HD=Hd-dD Union Hs-sD
  SZ =createDependencyMatrix(zc_jac) 
  dZ =createDependencyMatrix(ZC_discJac) # temp dependency to be used to determine HD2 and HZ2
  HZ1HD1=createDependencyToEventsDiscr(dD,dZ,evDep) 
  HZ2HD2=createDependencyToEventsCont(SD,SZ,evDep) 
  HZ=unionDependency(HZ1HD1[1],HZ2HD2[1])
  HD=unionDependency(HZ1HD1[2],HZ2HD2[2]) =#
  
 #=  zcf = Vector{Function}()
  for i = 1:length(odep.zceqs)# later change to numberZC
    push!(zcf, @RuntimeGeneratedFunction(odep.zceqs[i].args[2])) #args[2] cuz there is extra stuff
  end
  eventf = Vector{Function}()
  for i = 1:length(odep.eventEqus)# later change to numEvents
    push!(eventf, @RuntimeGeneratedFunction(odep.eventEqus[i].args[2])) 
  end  =# 
#######################################compute initial values##################################################
n=1
for k = 1:O # compute initial derivatives for x and q (similar to a recursive way )
  n=n*k
   for i = 1:T q[i].coeffs[k] = x[i].coeffs[k] end # q computed from x and it is going to be used in the next x
   for i = 1:T
      clearCache(taylorOpsCache,cacheSize);f(i,-1,-1,q,d, t ,taylorOpsCache)
      ndifferentiate!(integratorCache,taylorOpsCache[1] , k - 1)
      x[i].coeffs[k+1] = (integratorCache.coeffs[1]) / n # /fact cuz i will store der/fac like the convention...to extract the derivatives (at endof sim) multiply by fac  derderx=coef[3]*fac(2)
    end
end

for i = 1:T
  savedVars[i][1].coeffs .= x[i].coeffs  #to be changed  1 2 3 ?
  quantum[i] = relQ * abs(x[i].coeffs[1]) ;quantum[i]=quantum[i] < absQ ? absQ : quantum[i];quantum[i]=quantum[i] > maxErr ? maxErr : quantum[i] 
  computeNextTime(Val(O), i, initTime, nextStateTime, x, quantum)
  initSmallAdvance=0.1
  #t[0]=initTime#initSmallAdvance
  clearCache(taylorOpsCache,cacheSize);
  #@timeit "f" 
  f(i,-1,-1,q,d,t,taylorOpsCache);#@show taylorOpsCache
  computeNextInputTime(Val(O), i, initTime, initTime,taylorOpsCache[1] , nextInputTime, x,  quantum)
end

#@show nextStateTime,nextInputTime
for i=1:Z
  #= clearCache(taylorOpsCache,cacheSize);output=zcf[i](x,d,t,taylorOpsCache).coeffs[1]  =#
  clearCache(taylorOpsCache,cacheSize);
  #@timeit "zcf" 
  f(-1,i,-1,x,d,t,taylorOpsCache)        
                 
  oldsignValue[i,2]=taylorOpsCache[1][0] #value
  oldsignValue[i,1]=sign(taylorOpsCache[1][0]) #sign modify 
  #@timeit "compEvent" 
  computeNextEventTime(i,taylorOpsCache[1],oldsignValue,initTime,  nextEventTime, quantum)
  #@show output
end
#@show nextEventTime
#= @show x
@show q
@show quantum
@show nextStateTime =#
###################################################################################################################################################################
####################################################################################################################################################################
#---------------------------------------------------------------------------------while loop-------------------------------------------------------------------------
###################################################################################################################################################################
####################################################################################################################################################################
simt = initTime;count = 1 ;len=length(savedTimes);printcount=0;simulStepCount=0;totalSteps=0
  prevStepTime=initTime;prevStepVal = zeros(MVector{T,MVector{O+1,Float64}})
  for i = 1:T prevStepVal[i] .= x[i].coeffs end
  direction= zeros(MVector{T,Float64})
  flag= zeros(MVector{T,Float64})
  breakloop= zeros(MVector{1,Float64})
  maxIter=0
while simt < ft && maxIter < 50000
 # @show abort
 #= if breakloop
  break
 end =#
 if breakloop[1]>2.0
  break
end
  maxIter+=1
  sch = updateScheduler(nextStateTime,nextEventTime, nextInputTime)
  simt = sch[2]
  if  simt>ft  
   # println("sim ends at simt= ",simt)
    count += 1
    if len<count
      len=count*2
      for i=1:T
        resize!(savedVars[i],len)
        for z=count:len savedVars[i][z]=Taylor0(zeros(O+1),O) end# without this, the new ones are undefined
      end
      resize!(savedTimes,len)
    end
    for k = 1:T 
      integrateState(Val(O),x[k],integratorCache,ft-prevStepTime)
      savedVars[k][count].coeffs .=x[k].coeffs  
    end
    savedTimes[count]=ft#since simt passed ft, we could later interpolate instead
    break   ###################################################break##########################################
  end
  index = sch[1]
  numSteps[index]+=1
  t[0]=simt
  ##########################################state######################################## 
  if sch[3] == :ST_STATE
    #= if printcounter>0
      println("state step x= ",x)
      printcounter-=1
      #= if printcounter==0
        breakloop=true
      end =#
    end =#

   #=  if abs(q[index][0] - x[index][0]) > 10*quantum[index]
      breakloop[1]+=1.0
      @show index,simt,quantum[index]
      @show q[index][0],q[index][1]
      @show x[index][0],x[index][1],x[index][2]
      @show nextStateTime
  
    end =#


  
    elapsed = simt - tx[index];integrateState(Val(O),x[index],integratorCache,elapsed);tx[index] = simt 
    quantum[index] = relQ * abs(x[index].coeffs[1]) ;quantum[index]=quantum[index] < absQ ? absQ : quantum[index];quantum[index]=quantum[index] > maxErr ? maxErr : quantum[index]   
    for k = 1:O q[index].coeffs[k] = x[index].coeffs[k] end; tq[index] = simt    
    computeNextTime(Val(O), index, simt, nextStateTime, x, quantum) #
   # if simt>=3.09e-5
 
  #  end
    for i = 1:length(SD[index])
      j = SD[index][i] 
      if j != 0      
        elapsedx = simt - tx[j];if elapsedx > 0
         x[j].coeffs[1] = x[j](elapsedx);tx[j] = simt 
        # quantum[j] = relQ * abs(x[j].coeffs[1]) ;quantum[j]=quantum[j] < absQ ? absQ : quantum[j];quantum[j]=quantum[j] > maxErr ? maxErr : quantum[j]   
        end
        elapsedq = simt - tq[j];if elapsedq > 0 integrateState(Val(O-1),q[j],integratorCache,elapsedq);tq[j] = simt  end#q needs to be updated here for recomputeNext  
       
        #= if abs(q[j][0] - x[j][0]) > 10*quantum[j]
          breakloop[1]+=1.0
          @show j,simt,quantum[j]
          @show q[j][0],q[j][1]
          @show x[j][0],x[j][1],x[j][2]
          @show nextStateTime[j]
      
        end =#
                 
        for b in (jac[j]  )    
          elapsedq = simt - tq[b]
          if elapsedq>0
            integrateState(Val(O-1),q[b],integratorCache,elapsedq);tq[b]=simt
          end
       # end
      end
        clearCache(taylorOpsCache,cacheSize);f(j,-1,-1,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[j], taylorOpsCache[1],integratorCache,elapsed)
      
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
               #=  println("after recompute")
                @show nextStateTime =#
     #   end
       # println("$j under normal dependency")
      end#end if j!=0
    end#end for SD
   #=  for i = 1:length(SZ[index])
      j = SZ[index][i] 
      if j != 0 
        for b = 1:T # elapsed update all other vars that this derj depends upon.
          if zc_SimpleJac[j][b] != 0     
            elapsedx = simt - tx[b];if elapsedx>0 integrateState(Val(O),x[b],integratorCache,elapsedx);tx[b]=simt end
            #elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],integratorCache,elapsedq);tq[b]=simt end
          end
        end            
        #= clearCache(taylorOpsCache,cacheSize)#normally and later i should update x or q (integrate q=q+e derQ  for higher orders)
        computeNextEventTime(j,zcf[j](x,d,t,taylorOpsCache),oldsignValue,simt,  nextEventTime, quantum) =#
        clearCache(taylorOpsCache,cacheSize);f(-1,j,-1,x,d,t,taylorOpsCache)        
        computeNextEventTime(j,taylorOpsCache[1],oldsignValue,simt,  nextEventTime, quantum)
        #= if 0.4>simt > 0.31
        println("$index $j nexteventtime from SZ= ",nextEventTime)
        @show x[index]
        end =#
      end  #end if j!=0
    end#end for SZ
    ##################################input########################################
  elseif sch[3] == :ST_INPUT  # time of change has come to a state var that does not depend on anything...no one will give you a chance to change but yourself    
    #@show index,x
    if printcounter>0
      #println("input step index= ",index,"x= ",x)
      printcounter-=1
    end
    elapsed = simt - tx[index];integrateState(Val(O),x[index],integratorCache,elapsed);tx[index] = simt 
    quantum[index] = relQ * abs(x[index].coeffs[1]) ;quantum[index]=quantum[index] < absQ ? absQ : quantum[index];quantum[index]=quantum[index] > maxErr ? maxErr : quantum[index]   
    for k = 1:O q[index].coeffs[k] = x[index].coeffs[k] end; tq[index] = simt 
    clearCache(taylorOpsCache,cacheSize);f(index,-1,-1,q,d,t,taylorOpsCache)
    computeNextInputTime(Val(O), index, simt, elapsed,taylorOpsCache[1] , nextInputTime, x,  quantum)
    computeDerivative(Val(O), x[index], taylorOpsCache[1],integratorCache,elapsed)
   # reComputeNextTime(Val(O), index, simt, nextStateTime, x, q, quantum)
   #if 0.4>simt > 0.31
  # if printcounter>0 println("$index  -nextInputTime from input= ",nextInputTime);@show x;@show simt end
    # end 
    for i = 1:length(SD[index])
      j = SD[index][i] 
      if j != 0      
        elapsedx = simt - tx[j];
        if elapsedx > 0 
          x[j].coeffs[1] = x[j](elapsedx);tx[j] = simt 
          quantum[j] = relQ * abs(x[j].coeffs[1]) ;quantum[j]=quantum[j] < absQ ? absQ : quantum[j];quantum[j]=quantum[j] > maxErr ? maxErr : quantum[j]   
        end
        elapsedq = simt - tq[j];if elapsedq > 0 integrateState(Val(O-1),q[j],integratorCache,elapsedq);tq[j] = simt  end#q needs to be updated here for recomputeNext                 
        for b = 1:T # elapsed update all other vars that this derj depends upon.
          if jac[j][b] != 0     
            elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],integratorCache,elapsedq);tq[b]=simt end
          end
        end
        clearCache(taylorOpsCache,cacheSize);f(j,-1,-1,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[j], taylorOpsCache[1],integratorCache,elapsed)
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
      end#end if j!=0
    end#end for
    for i = 1:length(SZ[index])
      j = SZ[index][i] 
      if j != 0   
        for b = 1:T # elapsed update all other vars that this derj depends upon.
          if zc_SimpleJac[j][b] != 0     
            elapsedx = simt - tx[b];if elapsedx>0 integrateState(Val(O),x[b],integratorCache,elapsedx);tx[b]=simt end
           # elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],integratorCache,elapsedq);tq[b]=simt end
          end
        end              
       #=  clearCache(taylorOpsCache,cacheSize)#normally and later i should update x,q (integrate q=q+e derQ  for higher orders)
        computeNextEventTime(j,zcf[j](x,d,t,taylorOpsCache),oldsignValue,simt,  nextEventTime, quantum)#,maxIterer)  =#
        clearCache(taylorOpsCache,cacheSize);f(-1,j,-1,x,d,t,taylorOpsCache)        
                 computeNextEventTime(j,taylorOpsCache[1],oldsignValue,simt,  nextEventTime, quantum)
      end  
    end
  #################################################################event########################################
  else
          printcounter=5
         #=  println("x at start of event")
          @show x
          @show simt 
          @show quantum =#
      
          for b = 1:T # elapsed update all other vars that this zc depends upon.
            if zc_SimpleJac[index][b] != 0     
              elapsedx = simt - tx[b];
              if elapsedx>0 
                integrateState(Val(O),x[b],integratorCache,elapsedx);
               # println("should update var1 b=",b)
                tx[b]=simt 
              end
             elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],integratorCache,elapsedq);tq[b]=simt end
            end
          end    
          modifiedIndex=0#first we have a zc happened which corresponds to nexteventtime and index (one of zc) but we want also the sign in O to know ev+ or ev- 
         
          clearCache(taylorOpsCache,cacheSize);f(-1,index,-1,x,d,t,taylorOpsCache)     
         #=  println(" just after event")
         
          @show t 
          @show nextEventTime    =#
                
          if (taylorOpsCache[1][0])>=-1e-10       # sign is not needed here
            modifiedIndex=2*index-1   # the  event that just occured is at  this index
            #println("event1")
          else
            modifiedIndex=2*index
          end  
        # println("x before event= ",x)  
           
         for b=1:T
                if evDep[modifiedIndex].evContRHS[b]!==0.0  # use 3 signs or nan() function
                  elapsedx = simt - tx[b];
                  if elapsedx>0 
                    integrateState(Val(O),x[b],integratorCache,elapsedx);
                  #  println("var2 here b = ",b)
                    
                    tx[b]=simt 
                   end
                end
           end
           
         ###### eventf[modifiedIndex](x,d,t,taylorOpsCache) #if a choice to use x instead of q in events, then i think there should be a q update after the eventexecuted
        # @show x,modifiedIndex
         f(-1,-1,modifiedIndex,x,d,t,taylorOpsCache)
        
          for i=1:T
            #------------event influences a Continete var: x already updated in event...here update quantum and q and computenext
            if evDep[modifiedIndex].evCont[i]!==NaN   # use 3 signs or nan() function
                #quantum[i] = relQ * abs(x[i].coeffs[1]) ;quantum[i]=quantum[i] < absQ ? absQ : quantum[i];quantum[i]=quantum[i] > maxErr ? maxErr : quantum[i] 
                q[i][0]=x[i][0];tx[i] = simt;tq[i] = simt # for liqss updateQ?
              #  computeNextTime(Val(O), i, simt, nextStateTime, x, quantum) 
            end
          end
          nextEventTime[index]=Inf   #investigate more 
          for i = 1:length(HD[modifiedIndex]) # care about dependency to this event only
            j = HD[modifiedIndex][i] #  
            if j != 0      
              elapsedx = simt - tx[j];if elapsedx > 0 x[j].coeffs[1] = x[j](elapsedx);tx[j] = simt;#= @show j,x[j] =# end
              elapsedq = simt - tq[j];if elapsedq > 0 integrateState(Val(O-1),q[j],integratorCache,elapsedq);tq[j] = simt;#= @show q[j] =#  end#q needs to be updated here for recomputeNext                 
              for b = 1:T # elapsed update all other vars that this derj depends upon.
                if jac[j][b] != 0     
                  elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],integratorCache,elapsedq);tq[b]=simt;#= @show q[b] =# end
                end
              end
              clearCache(taylorOpsCache,cacheSize);f(j,-1,-1,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[j], taylorOpsCache[1],integratorCache,elapsedx)
              reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
              #@show j,x
            end#end if j!=0
          end
          for i = 1:length(HZ[modifiedIndex])
                j = HZ[modifiedIndex][i] 
                  if j != 0   
                    for b = 1:T # elapsed update all other vars that this derj depends upon.
                      if zc_SimpleJac[j][b] != 0     
                        elapsedx = simt - tx[b];if elapsedx>0 integrateState(Val(O),x[b],integratorCache,elapsedx);tx[b]=simt end
                       #elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],integratorCache,elapsedq);tq[b]=simt end
                      end
                    end            
                  #= clearCache(taylorOpsCache,cacheSize) #normally and later i should update q (integrate q=q+e derQ  for higher orders)          
                  computeNextEventTime(j,zcf[j](x,d,t,taylorOpsCache),oldsignValue,simt,  nextEventTime, quantum)#,maxIterer) =#
                  
                  clearCache(taylorOpsCache,cacheSize);f(-1,j,-1,x,d,t,taylorOpsCache)        
                 computeNextEventTime(j,taylorOpsCache[1],oldsignValue,simt,  nextEventTime, quantum)
                end  
              # if 0.4>simt > 0.31  println("$index $j nexteventtime from HZ= ",nextEventTime)   end   
          end
         #=  println("x end of step event")
         @show x 
         @show q =# =#
  end#end state/input/event
  if simt > savetime || sch[3] ==:ST_EVENT
    count += 1
    #savetime += savetimeincrement #next savetime
    if len<count
      len=count*2
      for i=1:T
        resize!(savedVars[i],len)
        for z=count:len
        savedVars[i][z]=Taylor0(zeros(O+1),O) # without this, the new ones are undefined
        end
      end
      resize!(savedTimes,len)
    end
    if savedTimes[count-1]!=prevStepTime  #if last point has not already been saved             
        savedTimes[count]=prevStepTime
        for k = 1:T             
          savedVars[k][count].coeffs .=prevStepVal[k]
        end
        count += 1
        if len<count
          len=count*2
          for i=1:T
            resize!(savedVars[i],len)
            for z=count:len
            savedVars[i][z]=Taylor0(zeros(O+1),O) # without this, the new ones are undefined
            end
          end
          resize!(savedTimes,len)
        end
    end
    for k = 1:T     
      elapsed = simt - tx[k];integrateState(Val(O),x[k],integratorCache,elapsed);tx[k] = simt #in case this point did not get updated.
        savedVars[k][count].coeffs .=x[k].coeffs 
    end
    savetime += savetimeincrement #next savetime
    savedTimes[count]=simt
  end#end if save
  prevStepTime=simt
  for k = 1:T    #store prev temporarily         
  prevStepVal[k] .=x[k].coeffs 
  end
end#end while
for i=1:T# throw away empty points
resize!(savedVars[i],count)
end
resize!(savedTimes,count)
Sol(O,savedTimes, savedVars,"qss$O",string(nameof(f)),numSteps,absQ,count,0)#0 I track simulSteps 
end#end integrate

