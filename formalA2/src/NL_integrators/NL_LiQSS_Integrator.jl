function LiQSS_integrate(::Val{O}, s::LiQSS_data{T,Z,O}, odep::NLODEProblem{T,D,Z,Y},f::Function) where {O,T,D,Z,Y}
    #*********************************settings*****************************************
    ft = s.finalTime;initTime = s.initialTime;relQ = s.dQrel;absQ = s.dQmin;maxErr=s.maxErr;savetimeincrement=s.savetimeincrement;savetime = savetimeincrement
    #*********************************qss method data*****************************************
    quantum = s.quantum;nextStateTime = s.nextStateTime;nextEventTime = s.nextEventTime;nextInputTime = s.nextInputTime
    tx = s.tx;tq = s.tq;x = s.x;q = s.q;t=s.t
    savedVars=s.savedVars;savedTimes=s.savedTimes;integratorCache=s.integratorCache;taylorOpsCache=s.taylorOpsCache;cacheSize=odep.cacheSize
    a=s.initJac;u=s.u;tu=s.tu
    #*********************************problem info*****************************************
    d = odep.discreteVars
    jacobian = odep.jacobian
    discJac = odep.discreteJacobian
    zc_jac = odep.ZC_jacobian
    ZC_discJac = odep.ZC_jacDiscrete
    evDep = odep.eventDependencies
    #********************************helper values*******************************  
    qaux=s.qaux
    olddx=s.olddx
    numSteps = zeros(MVector{T,Int})
    oldsignValue = MMatrix{Z,2}(zeros(Z*2))  #usedto track if zc changed sign; each zc has a value and a sign 
    jac=changeBasicToInts(jacobian)# change from type nonisbits to int so that access is cheaper down
    #*******************************create dependencies**************************
    SD = createDependencyMatrix(jacobian)
    dD =createDependencyMatrix(discJac) # temp dependency to be used to determine HD1 and HZ1 HD=Hd-dD Union Hs-sD
    SZ =createDependencyMatrix(zc_jac) 
    dZ =createDependencyMatrix(ZC_discJac) # temp dependency to be used to determine HD2 and HZ2
    HZ1HD1=createDependencyToEventsDiscr(dD,dZ,evDep) 
    HZ2HD2=createDependencyToEventsCont(SD,SZ,evDep) 
    HZ=unionDependency(HZ1HD1[1],HZ2HD2[1])
    HD=unionDependency(HZ1HD1[2],HZ2HD2[2])

    zcf = Vector{Function}()
    for i = 1:length(odep.zceqs)# later change to numberZC
      push!(zcf, @RuntimeGeneratedFunction(odep.zceqs[i].args[2])) #args[2] cuz there is extra stuff
    end
    eventf = Vector{Function}()
    for i = 1:length(odep.eventEqus)# later change to numEvents
      push!(eventf, @RuntimeGeneratedFunction(odep.eventEqus[i].args[2])) 
    end  
    #######################################compute initial values##################################################
    n=1
    for k = 1:O # compute initial derivatives for x and q (similar to a recursive way )
      n=n*k
       for i = 1:T
          q[i].coeffs[k] = x[i].coeffs[k]  # q computed from x and it is going to be used in the next x
        end
        for i = 1:T
          clearCache(taylorOpsCache,cacheSize);f(i,q,d, t ,taylorOpsCache)
          ndifferentiate!(integratorCache,taylorOpsCache[1] , k - 1)
          x[i].coeffs[k+1] = (integratorCache.coeffs[1]) / n # /fact cuz i will store der/fac like the convention...to extract the derivatives (at endof sim) multiply by fac  derderx=coef[3]*fac(2)
        end
    end
    for i = 1:T
      p=1    
      for k=1:O# deleting this causes scheduler error
        p=p*k
        m=p/k
          u[i][i][k]=p*x[i][k]-m*q[i][k-1]*a[i][i] #  later we will investigate inconsistencies of using data stored vs */ factorial!!! ...also do not confuse getindex for taylor...[0] first element and u[i][1]...first element    ########||||||||||||||||||||||||||||||||||||liqss|||||||||||||||||||||||||||||||||||||||||
      end
    end
    for i = 1:T    
      savedVars[i][1].coeffs .= x[i].coeffs  
      quantum[i] = relQ * abs(x[i].coeffs[1]) ;quantum[i]=quantum[i] < absQ ? absQ : quantum[i];quantum[i]=quantum[i] > maxErr ? maxErr : quantum[i] 
     updateQ(Val(O),i,x,q,quantum,a,u,qaux,olddx,tx,tq,tu,initTime,ft,nextStateTime) 
    end
    for i = 1:T 
      clearCache(taylorOpsCache,cacheSize);f(i,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[i], taylorOpsCache[1],integratorCache,0.0)
      Liqss_reComputeNextTime(Val(O), i, initTime, nextStateTime, x, q, quantum,a)
      computeNextInputTime(Val(O), i, initTime, 0.1,taylorOpsCache[1] , nextInputTime, x,  quantum)        #   *******************************testing: remove comments later************************
    end
    for i=1:Z
      clearCache(taylorOpsCache,cacheSize)
      output=zcf[i](x,d,t,taylorOpsCache).coeffs[1] #test this evaluation
      oldsignValue[i,2]=output #value
      oldsignValue[i,1]=sign(output) #sign modify 
      computeNextEventTime(i,output,oldsignValue,initTime,  nextEventTime, quantum)#,printCounter)
    end
    ###################################################################################################################################################################
    ####################################################################################################################################################################
    #---------------------------------------------------------------------------------while loop-------------------------------------------------------------------------
    ###################################################################################################################################################################
    ####################################################################################################################################################################
simt = initTime;count = 1 ;len=length(savedTimes);printcount=0
prevStepTime=initTime;prevStepVal = zeros(MVector{T,MVector{O+1,Float64}})
for i = 1:T prevStepVal[i] .= x[i].coeffs end
while simt < ft && printcount < 500000
  printcount+=1
  sch = updateScheduler(nextStateTime,nextEventTime, nextInputTime)
  simt = sch[2]
  if  simt>ft  
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
        elapsed = simt - tx[index];integrateState(Val(O),x[index],integratorCache,elapsed);tx[index] = simt 
        quantum[index] = relQ * abs(x[index].coeffs[1]) ;quantum[index]=quantum[index] < absQ ? absQ : quantum[index];quantum[index]=quantum[index] > maxErr ? maxErr : quantum[index]    
        updateQ(Val(O),index,x,q,quantum,a,u,qaux,olddx,tx,tq,tu,simt,ft,nextStateTime) ;tq[index]=simt
        for i = 1:length(SD[index])
          j = SD[index][i] 
          if j != 0      
            elapsedx = simt - tx[j];if elapsedx > 0 x[j].coeffs[1] = x[j](elapsedx);tx[j] = simt end
            elapsedq = simt - tq[j];if elapsedq > 0 integrateState(Val(O-1),q[j],integratorCache,elapsedq);tq[j] = simt  end#q needs to be updated here for recomputeNext                 
            for b = 1:T # elapsed update all other vars that this derj depends upon.
              if jac[j][b] != 0     
                elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],integratorCache,elapsedq);tq[b]=simt end
              end
            end
            clearCache(taylorOpsCache,cacheSize);f(j,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[j], taylorOpsCache[1],integratorCache,elapsed)
            reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
          end#end if j!=0
        end#end for SD
        for i = 1:length(SZ[index])
          j = SZ[index][i] 
          if j != 0             
            #normally and later i should update q (integrate q=q+e derQ  for higher orders)
            clearCache(taylorOpsCache,cacheSize)
            computeNextEventTime(j,zcf[j](x,d,t,taylorOpsCache)[0],oldsignValue,simt,  nextEventTime, quantum)#,printCounter)
          end  #end if j!=0
        end#end for SZ
        updateLinearApprox(Val(O),index,x,q,a,u,qaux,olddx,tu,simt)########||||||||||||||||||||||||||||||||||||liqss|||||||||||||||||||||||||||||||||||||||||
        ##################################input########################################
  elseif sch[3] == :ST_INPUT  # time of change has come to a state var that does not depend on anything...no one will give you a chance to change but yourself    
        clearCache(taylorOpsCache,cacheSize);f(index,q,d,t,taylorOpsCache)
        elapsed = simt - tx[index];computeNextInputTime(Val(O), index, simt, elapsed,taylorOpsCache[1] , nextInputTime, x,  quantum)
        for i = 1:length(SD[index])
          j = SD[index][i] 
          if j != 0      
            elapsedx = simt - tx[j];if elapsedx > 0 x[j].coeffs[1] = x[j](elapsedx);tx[j] = simt end
            elapsedq = simt - tq[j];if elapsedq > 0 integrateState(Val(O-1),q[j],integratorCache,elapsedq);tq[j] = simt  end#q needs to be updated here for recomputeNext                 
            for b = 1:T # elapsed update all other vars that this derj depends upon.
              if jac[j][b] != 0     
                elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],integratorCache,elapsedq);tq[b]=simt end
              end
            end
            clearCache(taylorOpsCache,cacheSize);f(j,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[j], taylorOpsCache[1],integratorCache,elapsed)
            reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
          end#end if j!=0
        end#end for
        for i = 1:length(SZ[index])
          j = SZ[index][i] 
          if j != 0             
            clearCache(taylorOpsCache,cacheSize)#normally and later i should update x,q (integrate q=q+e derQ  for higher orders)
            computeNextEventTime(j,zcf[j](q,d,t,taylorOpsCache)[0],oldsignValue,simt,  nextEventTime, quantum)#,printCounter) 
          end  
        end
      #################################################################event########################################
  else
        modifiedIndex=0#first we have a zc happened which corresponds to nexteventtime and index (one of zc) but we want also the sign in O to know ev+ or ev- 
        if (zcf[index](x,d,t,taylorOpsCache).coeffs[1])>0       # sign is not needed here
          modifiedIndex=2*index-1   # the  event that just occured is at  this index
        else
          modifiedIndex=2*index
        end       
        eventf[modifiedIndex](q,d,t,taylorOpsCache) #if a choice to use x instead of q in events, then i think there should be a q update after the eventexecuted
        nextEventTime[index]=Inf   #investigate more
        for i = 1:length(HD[modifiedIndex]) # care about dependency to this event only
          j = HD[modifiedIndex][i] #  
          if j != 0      
            elapsedx = simt - tx[j];if elapsedx > 0 x[j].coeffs[1] = x[j](elapsedx);tx[j] = simt end
            elapsedq = simt - tq[j];if elapsedq > 0 integrateState(Val(O-1),q[j],integratorCache,elapsedq);tq[j] = simt  end#q needs to be updated here for recomputeNext                 
            for b = 1:T # elapsed update all other vars that this derj depends upon.
              if jac[j][b] != 0     
                elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],integratorCache,elapsedq);tq[b]=simt end
              end
            end
            clearCache(taylorOpsCache,cacheSize);f(j,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[j], taylorOpsCache[1],integratorCache,elapsed)
            reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
          end#end if j!=0
        end
        for i = 1:length(HZ[modifiedIndex])
              j = HZ[modifiedIndex][i] 
                if j != 0             
                clearCache(taylorOpsCache,cacheSize) #normally and later i should update q (integrate q=q+e derQ  for higher orders)          
                computeNextEventTime(j,zcf[j](x,d,t,taylorOpsCache)[0],oldsignValue,simt,  nextEventTime, quantum)#,printCounter)
              end        
        end
  end#end state/input/event
  if simt > savetime
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
      if savedTimes[count-1]!=prevStepTime  #if last point has not already been saved             
          savedTimes[count]=prevStepTime
          for k = 1:T             
            savedVars[k][count].coeffs .=prevStepVal[k]
          end
          count += 1
      end
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
Sol(O,savedTimes, savedVars,"liqss$O",string(nameof(f)),numSteps,absQ,0)
end#end integrate
    
    