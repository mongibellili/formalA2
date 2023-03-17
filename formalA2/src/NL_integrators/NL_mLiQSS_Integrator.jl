 #using TimerOutputs
 function mLiQSS_integrate(::Val{O}, s::LiQSS_data{T,Z,O}, odep::NLODEProblem{T,D,Z,Y},f::Function) where {O,T,D,Z,Y}
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
qaux=s.qaux;olddx=s.olddx;olddxSpec = zeros(MVector{T,MVector{O,Float64}}) # later can only care about 1st der
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
  p=1
  for k=1:O
    p=p*k
    m=p/k
    for j=1:T
      if j!=i
        u[i][j][k]=p*x[i][k]-a[i][i]*m*q[i][k-1]-a[i][j]*m*q[j][k-1]
      else
        u[i][j][k]=p*x[i][k]-a[i][i]*m*q[i][k-1]
      end
    end
  end
end
for i = 1:T
  savedVars[i][1].coeffs .= x[i].coeffs  #to be changed  1 2 3 ?
  quantum[i] = relQ * abs(x[i].coeffs[1]) ;quantum[i]=quantum[i] < absQ ? absQ : quantum[i];quantum[i]=quantum[i] > maxErr ? maxErr : quantum[i] 
  updateQ(Val(O),i,x,q,quantum,a,u,qaux,olddx,tx,tq,tu,initTime,ft,nextStateTime) 
end
for i = 1:T
  clearCache(taylorOpsCache,cacheSize);f(i,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[i], taylorOpsCache[1],integratorCache,0.0)#0.0 used to be elapsed...even down below not neeeded anymore
  Liqss_reComputeNextTime(Val(O), i, initTime, nextStateTime, x, q, quantum,a)
  computeNextInputTime(Val(O), i, initTime, 0.1,taylorOpsCache[1] , nextInputTime, x,  quantum)#not complete, currently elapsed=0.1 is temp until fixed
end
for i=1:Z
  clearCache(taylorOpsCache,cacheSize)
  output=zcf[i](x,d,t,taylorOpsCache).coeffs[1] #test this evaluation
  oldsignValue[i,2]=output #value
  oldsignValue[i,1]=sign(output) #sign modify 
  computeNextEventTime(i,output,oldsignValue,initTime,  nextEventTime, quantum)#,printCounter)# if Z=0 nexteventtime did not get filled with zeros in the first place
end
###################################################################################################################################################################
####################################################################################################################################################################
#---------------------------------------------------------------------------------while loop-------------------------------------------------------------------------
###################################################################################################################################################################
#################################################################################################################################################################### 
simt = initTime;count = 1 ;len=length(savedTimes);printcount=0;simulStepCount=0
prevStepTime=initTime;prevStepVal = zeros(MVector{T,MVector{O+1,Float64}})
for i = 1:T prevStepVal[i] .= x[i].coeffs end
while simt < ft && printcount < 80000000
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
      updateQ(Val(O),index,x,q,quantum,a,u,qaux,olddx,tx,tq,tu,simt,ft,nextStateTime) ;tq[index] = simt        
      #----------------------------------------------------check dependecy cycles---------------------------------------------      
      for l = 1:length(SD[index])
      j = SD[index][l] 
      if j != 0 && j!=index && a[index][j]*a[j][index]!=0                  
        olddxSpec[index][1]=x[index][1]
        if isCycle_and_simulUpdate(Val(O),index,j,x,q,quantum,a,u,qaux,olddx,olddxSpec,tx,tq,tu,simt,ft)
          simulStepCount+=1             
          for b = 1:T # elapsed update all other vars that these der i & j depend upon.needed for when sys has 3 or more vars.
            if jac[j][b] != 0  ||  jac[index][b] != 0      
              elapsedq = simt - tq[b] ;if elapsedq>0 integrateState(Val(O-1),q[b],integratorCache,elapsedq);tq[b]=simt end
            end
          end
          #compute olddxSpec_i using new qi and qjaux to annihilate the influence of qi (keep j influence only) when finding aij=(dxi-dxi)/(qj-qjaux)
          qjtemp=q[j][0];q[j][0]=qaux[j][1]
          clearCache(taylorOpsCache,cacheSize);f(index,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[index], taylorOpsCache[1],integratorCache,elapsed)
          clearCache(taylorOpsCache,cacheSize);f(j,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[j], taylorOpsCache[1],integratorCache,elapsed)
          olddx[j][1]=x[j][1]  # needed to find a_jj (qi annihilated, qj kept)
          olddxSpec[index][1]= x[index][1] # new qi used now so it does not have an effect later on aij
          q[j][0]=qjtemp  # get back qj

          qitemp=q[index][0];q[index][0]=qaux[index][1]# 
          clearCache(taylorOpsCache,cacheSize);f(index,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[index], taylorOpsCache[1],integratorCache,elapsed)
          clearCache(taylorOpsCache,cacheSize);f(j,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[j], taylorOpsCache[1],integratorCache,elapsed)
          olddx[index][1]=x[index][1]              
          olddxSpec[j][1]=x[j][1]  # new qj used now so it does not have an effect later on aji
          q[index][0]=qitemp

          clearCache(taylorOpsCache,cacheSize);f(index,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[index], taylorOpsCache[1],integratorCache,elapsed)
          clearCache(taylorOpsCache,cacheSize);f(j,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[j], taylorOpsCache[1],integratorCache,elapsed)
          Liqss_reComputeNextTime(Val(O), index, simt, nextStateTime, x, q, quantum,a)
          Liqss_reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum,a)

          updateOtherApprox(Val(O),index,j,x,q,a,u,qaux,olddxSpec,tu,simt)
          updateOtherApprox(Val(O),j,index,x,q,a,u,qaux,olddxSpec,tu,simt)

          for l = 1:length(SD[j])
                k = SD[j][l] 
                if k != 0  && k!=index && k!=j
                  elapsedx = simt - tx[k]
                  if elapsedx > 0
                    x[k].coeffs[1] = x[k](elapsedx);tx[k] = simt
                    differentiate!(integratorCache,x[k])
                    x[k][1] = integratorCache(elapsedx)
                    olddxSpec[k][1]=x[j][1]
                  end
                  elapsedq = simt - tq[k]
                  if elapsedq > 0     
                    integrateState(Val(O-1),q[k],integratorCache,elapsedq); tq[k] = simt
                  end
                  for b = 1:T # elapsed update all other vars that this derj depends upon.needed for when sys has 3 or more vars.
                    if jac[k][b] != 0       
                      elapsedq = simt - tq[b]
                      if elapsedq>0
                        integrateState(Val(O-1),q[b],integratorCache,elapsedq);tq[b]=simt
                      end
                    end
                  end                     
                  clearCache(taylorOpsCache,cacheSize);f(k,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[k], taylorOpsCache[1],integratorCache,elapsed)
                  Liqss_reComputeNextTime(Val(O), k, simt, nextStateTime, x, q, quantum,a)
                  updateOtherApprox(Val(O),k,j,x,q,a,u,qaux,olddxSpec,tu,simt)
                end#end if k!=0
          end#end for k depend on j
          for l = 1:length(SZ[j])
                k = SZ[index][l] 
                if k != 0             
                  #normally and later i should update q (integrate q=q+e derQ  for higher orders)
                  clearCache(taylorOpsCache,cacheSize)
                  computeNextEventTime(k,zcf[k](x,d,t,taylorOpsCache)[0],oldsignValue,simt,  nextEventTime, quantum)#,printCounter)
                end  #end if j!=0
          end#end for SZ                                        
          updateLinearApprox(Val(O),j,x,q,a,u,qaux,olddx,tu,simt)             
        end#end ifcycle check
        # tx[j] = simt #  for sys with 2 vars: i won't be updated because of elapsed and j wont be updated
        # tq[j] = simt
      end
    end#end FOR_cycle check
    #-------------------------------------------------------------------------------------
    #---------------------------------normal liqss: proceed--------------------------------
    #-------------------------------------------------------------------------------------

    for l = 1:length(SD[index])
      j = SD[index][l] 
      if j != 0           
        elapsedx = simt - tx[j]
        if elapsedx > 0
          x[j].coeffs[1] = x[j](elapsedx);tx[j] = simt
          differentiate!(integratorCache,x[j])
          x[j][1] = integratorCache(elapsedx)
          olddxSpec[j][1]=x[j][1] # if elapsedx>0 then elapsedq>0 (confirm?)
        end
        elapsedq = simt - tq[j]
        if elapsedq > 0
          integrateState(Val(O-1),q[j],integratorCache,elapsedq);tq[j] = simt         
        end
        for b = 1:T # elapsed update all other vars that this derj depends upon
          if jac[j][b] != 0      
            elapsedq = simt - tq[b]
            if elapsedq>0
              integrateState(Val(O-1),q[b],integratorCache,elapsedq);tq[b]=simt
            end
          end
        end
        clearCache(taylorOpsCache,cacheSize);f(j,q,d,t,taylorOpsCache)
        if x[j][1]!=taylorOpsCache[1][0]#if none of the above q changed then the der would be same and no need for wasting resources
          computeDerivative(Val(O), x[j], taylorOpsCache[1],integratorCache,elapsed)
          Liqss_reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum,a)
          updateOtherApprox(Val(O),j,index,x,q,a,u,qaux,olddxSpec,tu,simt)# Later inverstigate oldspec not updated when aij*aji=0 above
        end
      end#end if j!=0
    end#end for SD
    for l = 1:length(SZ[index])
      j = SZ[index][l] 
      if j != 0             
        #normally and later i should update q (integrate q=q+e derQ  for higher orders)
        clearCache(taylorOpsCache,cacheSize)
        computeNextEventTime(j,zcf[j](x,d,t,taylorOpsCache)[0],oldsignValue,simt,  nextEventTime, quantum)#,printCounter)
      end  #end if j!=0
    end#end for SZ
    # if abs(a[index][index])>1e-6  # if index depends on itself update, otherwise leave zero 
    updateLinearApprox(Val(O),index,x,q,a,u,qaux,olddx,tu,simt)########||||||||||||||||||||||||||||||||||||liqss|||||||||||||||||||||||||||||||||||||||||
  #  end

    ##################################input########################################
  elseif sch[3] == :ST_INPUT  # time of change has come to a state var that does not depend on anything...no one will give you a chance to change but yourself  
    println("input event")
    elapsed = simt - tx[index]    
    clearCache(taylorOpsCache,cacheSize)   
    f(index,q,d,t,taylorOpsCache)
      computeNextInputTime(Val(O), index, simt, elapsed,taylorOpsCache[1] , nextInputTime, x,  quantum)
    for i = 1:length(SD[index])
      j = SD[index][i] 
      if j != 0             
        elapsed = simt - tx[j]
        if elapsed > 0

          x[j].coeffs[1] = x[j](elapsed)#.coeffs[1] #evaluate x at new time only...derivatives get updated next using computeDerivati()
          tx[j] = simt
        end
        quantum[j] = relQ * abs(x[j].coeffs[1]) #derx=coef[2]*fac(1), derderx=coef[3]*fac(2)            
        quantum[j]=quantum[j] < absQ ? absQ : quantum[j]
        quantum[j]=quantum[j] > maxErr ? maxErr : quantum[j] 
        clearCache(taylorOpsCache,cacheSize)
        f(j,q,d,t,taylorOpsCache)
        computeDerivative(Val(O), x[j], taylorOpsCache[1],integratorCache,elapsed)
        # computeDerivative(Val(O), x[j], taylorOpsCache[1])
        Liqss_reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum,a)
      end#end if j!=0
    end#end for
    for i = 1:length(SZ[index])
      j = SZ[index][i] 
      if j != 0             
        #normally and later i should update q (integrate q=q+e derQ  for higher orders)
        clearCache(taylorOpsCache,cacheSize)
        
        computeNextEventTime(j,zcf[j](q,d,t,taylorOpsCache)[0],oldsignValue,simt,  nextEventTime, quantum)#,printCounter) 
      end  
      # println("end input:who is resizing?")
    end
  #################################################################event########################################
  else
    #first we have a zc happened which corresponds to nexteventtime and index (one of zc) but we want also the sign in O to know ev+ or ev- 

    modifiedIndex=0
    if (zcf[index](x,d,t,taylorOpsCache).coeffs[1])>0       # sign is not needed here
      modifiedIndex=2*index-1   # the  event that just occured is at  this index
    else
      modifiedIndex=2*index
    end       
    eventf[modifiedIndex](q,d,t,taylorOpsCache)
    #if a choice to use x instead of q in events, then i think there should be a q update after the eventexecuted
    nextEventTime[index]=Inf   #investigate more
    for i = 1:length(HD[modifiedIndex]) # care about dependency to this event only
        j = HD[modifiedIndex][i] # 
        if j != 0
              elapsed = simt - tx[j]             
              if elapsed > 0  # if event triggere by change of sign and time=now then elapsed=0
                  # println("elapsed appear?= ",elapsed) #only appear when zc depends on x1 and it affects derx2 && derx1 doesn't depend on x1             
                  x[j].coeffs[1] = x[j](elapsed)#.coeffs[1] #evaluate x at new time only...derivatives get updated next using computeDerivativ()
                  q[j].coeffs[1] = x[j].coeffs[1]
                  tx[j] = simt
              end
                clearCache(taylorOpsCache,cacheSize)

              f(j,q,d,t,taylorOpsCache)
              computeDerivative(Val(O), x[j], taylorOpsCache[1],integratorCache,elapsed)
              # computeDerivative(Val(O), x[j], taylorOpsCache[1])
              Liqss_reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum,a)  
        end     
    end
    for i = 1:length(HZ[modifiedIndex])
          j = HZ[modifiedIndex][i] # this line before the next line or vice versa gave the same bench results
            if j != 0             
            #normally and later i should update q (integrate q=q+e derQ  for higher orders)
            clearCache(taylorOpsCache,cacheSize)           
            computeNextEventTime(j,zcf[j](x,d,t,taylorOpsCache)[0],oldsignValue,simt,  nextEventTime, quantum)#,printCounter)
          end        
    end
  end#end state/input/event

 
  
  if simt > savetime
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
        #if simt-savetime>2*savetimeincrement# in case there is a jump/gap, we want to save the last/prevoius pt (in addition )
          if savedTimes[count-1]!=prevStepTime  #if last point has not already been saved             
            savedTimes[count]=prevStepTime
            for k = 1:T             
              savedVars[k][count].coeffs .=prevStepVal[k]
            end
            count += 1
          end
        # end
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
  for k = 1:T             
    prevStepVal[k] .=x[k].coeffs 
  end

end#end while
for i=1:T# throw away empty points
  resize!(savedVars[i],count)
end
resize!(savedTimes,count)
Sol(O,savedTimes, savedVars,"mliqss$O",string(nameof(f)),numSteps,absQ,simulStepCount)
end#end integrate
    
    