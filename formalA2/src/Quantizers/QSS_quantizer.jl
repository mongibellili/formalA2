@inline function integrateState(::Val{0}, x::Taylor0{Float64},cacheT::Taylor0{Float64},elapsed::Float64) #@inline increased performance
  #nothing: created for elapse-updating q in order1 which does not happen
end

@inline function integrateState(::Val{1}, x::Taylor0{Float64},cacheT::Taylor0{Float64},elapsed::Float64) #@inline increased performance
  x.coeffs[1] = x(elapsed)
end
@inline function integrateState(::Val{2}, x::Taylor0{Float64},cacheT::Taylor0{Float64},elapsed::Float64) #@inline increased performance
 # @show x
  x.coeffs[1] = x(elapsed)
  #@show x
  differentiate!(cacheT,x)
  x.coeffs[2] = cacheT(elapsed)
  #@show x
  #cacheT.coeffs.=0.0 #clear cache
end
@inline function integrateState(::Val{3}, x::Taylor0{Float64},cacheT::Taylor0{Float64},elapsed::Float64) #@inline increased performance
  #@show elapsed
  #@show x
  x.coeffs[1] = x(elapsed)
  differentiate!(cacheT,x)
  
  x.coeffs[2] = cacheT(elapsed)
  ndifferentiate!(cacheT,x,2)
  x.coeffs[3] = cacheT(elapsed)/2
 
 #=  e=elapsed
    x[0]=x[0]+e*x[1]+e*e*2*x[2]/2+e*e*e*6*x[3]/6
  #  @show x[0]
    x[1]=x[1]+e*2*x[2]+e*e*6*x[3]/2
   # @show x[1]
    x[2]=(2*x[2]+e*6*x[3])/2 =#
 #@show x
end

######################################################################################################################################"
function computeDerivative( ::Val{1}  ,x::Taylor0{Float64},f::Taylor0{Float64},cache::Taylor0{Float64},elap::Float64  )
   #x.coeffs[2] =f(elap)
   x.coeffs[2] =f.coeffs[1]
   return nothing
end
function computeDerivative( ::Val{2}  ,x::Taylor0{Float64},f::Taylor0{Float64},cache::Taylor0{Float64},elap::Float64 ) # approx deriv ...some alloc gain in 11us for FT=5.0  58---69
   x.coeffs[2] =f.coeffs[1]
   x.coeffs[3] =f.coeffs[2]/2
   return nothing
end
#= function computeDerivative( ::Val{2}  ,x::Taylor0{Float64},f::Taylor0{Float64},cache::Taylor0{Float64},elap::Float64 )
  x.coeffs[2] =f(elap)
  differentiate!(cache,f)
  x.coeffs[3]=cache(elap)/2
  return nothing
end =#
#= function computeDerivative( ::Val{3}  ,x::Taylor0{Float64},f::Taylor0{Float64},cache::Taylor0{Float64},elap::Float64 )
  x.coeffs[2] =f(elap)
  differentiate!(cache,f)
  x.coeffs[3]=cache(elap)/2
  ndifferentiate!(cache,f,2)
  x.coeffs[4]=cache(elap)/6
  return nothing
end =#

function computeDerivative( ::Val{3}  ,x::Taylor0{Float64},f::Taylor0{Float64},cache::Taylor0{Float64},elap::Float64 )
  x.coeffs[2] =f[0]#f(elap)
  x.coeffs[3]=f.coeffs[2]/2
  x.coeffs[4]=f.coeffs[3]/3 # coeff3*2/6
  return nothing
end
######################################################################################################################################"
function computeNextTime(::Val{1}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, x::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T}#i can be absorbed
  absDeltaT=1e-12 # minimum deltaT to protect against der=Inf coming from sqrt(0) for example...similar to min ΔQ
    if (x[i].coeffs[2]) != 0
        tempTime=max(abs(quantum[i] /(x[i].coeffs[2])),absDeltaT)# i can avoid the use of max
        if tempTime!=absDeltaT #normal
            nextTime[i] = currentTime + tempTime#sqrt(abs(quantum[i] / ((x[i].coeffs[3])*2))) #*2 cuz coeff contains fact()
        else#usual (quant/der) is very small
          x[i].coeffs[2]=sign(x[i].coeffs[2])*(abs(quantum[i])/absDeltaT)# adjust  derivative if it is too high
          nextTime[i] = currentTime + tempTime
          #println("smalldelta")
        end
    else
      nextTime[i] = Inf
    end
    return nothing
end

function computeNextTime(::Val{2}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, x::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T}
    absDeltaT=1e-12 # minimum deltaT to protect against der=Inf coming from sqrt(0) for example...similar to min ΔQ
      if (x[i].coeffs[3]) != 0
          tempTime=max(sqrt(abs(quantum[i] / ((x[i].coeffs[3])))),absDeltaT)
          
          if tempTime!=absDeltaT #normal
              nextTime[i] = currentTime + tempTime#sqrt(abs(quantum[i] / ((x[i].coeffs[3])*2))) #*2 cuz coeff contains fact()
          else#usual sqrt(quant/der) is very small
            x[i].coeffs[3]=sign(x[i].coeffs[3])*(abs(quantum[i])/(absDeltaT*absDeltaT))/2# adjust second derivative if it is too high
            nextTime[i] = currentTime + tempTime
           # println("smalldelta in compute next")
          end
      else
        nextTime[i] = Inf
      end
      return nothing
end
function computeNextTime(::Val{3}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, x::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T}
  absDeltaT=1e-12 # minimum deltaT to protect against der=Inf coming from sqrt(0) for example...similar to min ΔQ
    if (x[i].coeffs[4]) != 0
        tempTime=max(cbrt(abs(quantum[i] / ((x[i].coeffs[4])))),absDeltaT)   #6/6
        if tempTime!=absDeltaT #normal
            nextTime[i] = currentTime + tempTime#sqrt(abs(quantum[i] / ((x[i].coeffs[3])*2))) #*2 cuz coeff contains fact()
        else#usual sqrt(quant/der) is very small
          x[i].coeffs[4]=sign(x[i].coeffs[4])*(abs(quantum[i])/(absDeltaT*absDeltaT*absDeltaT))/6# adjust third derivative if it is too high
          nextTime[i] = currentTime + tempTime
        end
    else
      nextTime[i] = Inf
    end
    return nothing
end
######################################################################################################################################"
function reComputeNextTime(::Val{1}, index::Int, currentTime::Float64, nextTime::MVector{T,Float64}, x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T}
  coef=@SVector [q[index].coeffs[1] - (x[index].coeffs[1]) - quantum[index], -x[index].coeffs[2]]
  time1 = currentTime + minPosRoot(coef, Val(1))
  coef=setindex(coef,q[index].coeffs[1] - (x[index].coeffs[1]) + quantum[index],1)
  time2 = currentTime + minPosRoot(coef, Val(1))
  nextTime[index] = time1 < time2 ? time1 : time2
  #later guard against very small Δt like in computenext
  return nothing
end

function reComputeNextTime(::Val{2}, index::Int, currentTime::Float64, nextTime::MVector{T,Float64}, x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T}
  absDeltaT=1e-12
  coef=@SVector [q[index].coeffs[1] - (x[index].coeffs[1]) - quantum[index], q[index].coeffs[2]-x[index].coeffs[2],-(x[index].coeffs[3])]#not *2 because i am solving c+bt+a/2*t^2
  time1 =  minPosRoot(coef, Val(2))
  coef=setindex(coef,q[index].coeffs[1] - (x[index].coeffs[1]) + quantum[index],1)
  time2 =  minPosRoot(coef, Val(2))
  timeTemp = time1 < time2 ? time1 : time2
  tempTime=max(timeTemp,absDeltaT)

  if tempTime==absDeltaT #normal
   # x[index].coeffs[3]=sign(x[index].coeffs[3])*(abs(quantum[index])/(absDeltaT*absDeltaT))/2
    println("smalldelta in recompute next")
  end



  nextTime[index] = currentTime +tempTime

  #later guard against very small Δt like in computenext
  return nothing
end
  
function reComputeNextTime(::Val{3}, index::Int, currentTime::Float64, nextTime::MVector{T,Float64}, x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T}
  #coef=@SVector [q[index].coeffs[1] - (x[index].coeffs[1]) - quantum[index], q[index].coeffs[2]-x[index].coeffs[2],(q[index].coeffs[3])-(x[index].coeffs[3]),-(x[index].coeffs[4])]
  #time1 = currentTime + minPosRoot(coef, Val(3))
  #pp=pointer(Vector{NTuple{2,Float64}}(undef, 5))
  a=-x[index][3];b=q[index][2]-x[index][2];c=q[index][1]-x[index][1];d=q[index][0]-x[index][0]-quantum[index]
  #GC.enable(false)
 #time1 = currentTime+smallestpositiverootintervalnewtonregulafalsi((a, b, c, d))
  time1 = currentTime+cubic5(a, b, c, d)
 # GC.enable(true)
  #coef=setindex(coef,q[index].coeffs[1] - (x[index].coeffs[1]) + quantum[index],1)
  d=q[index][0]-x[index][0]+quantum[index]
 # GC.enable(false)
  #time2 = currentTime+smallestpositiverootintervalnewtonregulafalsi((a, b, c, d))
  time2 = currentTime+cubic5(a, b, c, d)
 # GC.enable(true)
  #time2 = currentTime + minPosRoot(coef, Val(3))
  nextTime[index] = time1 < time2 ? time1 : time2
 #=  if x[index][1]==0.46926407070799725 && q[index][1]==0.4289402307079971
    @show a,b,c,d
    @show time1,time2
    @show nextTime[index]
  end =#
 #=  if nextTime[index]==Inf && x[index][1]!=0.0 && x[index][2]!=0.0 && x[index][3]!=0.0  # in oregonator i have x diverge from q cuz of bad root finders
    nextTime[index]=currentTime+1e-9
    #nextTime[index]=cbrt(abs(quantum[index] / x[index][3]))
    #nextTime[index]=sqrt(abs(quantum[index] / x[index][2]))
  end =#
  #later guard against very small Δt like in computenext
  return nothing
end
######################################################################################################################################"
function computeNextInputTime(::Val{1}, i::Int, currentTime::Float64,elapsed::Float64, tt::Taylor0{Float64} ,nextInputTime::MVector{T,Float64}, x::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T}
    df=0.0
    oldDerX=x[i].coeffs[2]
    newDerX=tt.coeffs[1] 
    if elapsed > 0
      df=(newDerX-oldDerX)/elapsed
    else
      df= quantum[i]*1e6
    end     
    if df!=0.0
       nextInputTime[i]=currentTime+sqrt(abs(2*quantum[i] / df))
    else
      nextInputTime[i] = Inf
    end
    return nothing
end

  
function computeNextInputTime(::Val{2}, i::Int, currentTime::Float64,elapsed::Float64, tt::Taylor0{Float64} ,nextInputTime::MVector{T,Float64}, x::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T}
  df=0.0
  oldDerDerX=((x[i].coeffs[3])*2.0)
  #@show x
  newDerDerX=(tt.coeffs[2])# 1st der of tt cuz tt itself is derx=f
  #@show newDerDerX
  if elapsed > 0.0
      df=(newDerDerX-oldDerDerX)/(elapsed)
     # println("df=new-old= ",df)
  else
      df= quantum[i]*1e6#*1e12
  end       
  if df!=0.0
      nextInputTime[i]=currentTime+cbrt((abs(quantum[i]/df)))    #df mimics 3rd der 
     # println("usedcbrt")
  else
      if newDerDerX==0 && x[i][1]!=0 #predicted second derivative is 0 & should be not be used to determine nexttime. use 1st der
         #nextInputTime[i]=currentTime+(abs(2*quantum[i] / x[i][1]))  #*2 is not analytic: is just there to increase stepsize
         nextInputTime[i]=currentTime+sqrt(abs(1*quantum[i] / x[i][1]))  #I used the same formulae even with 1st der so that it is fair to other vars
      else
          nextInputTime[i] = Inf
      end
  end
    return nothing
end
function computeNextInputTime(::Val{3}, i::Int, currentTime::Float64,elapsed::Float64, tt::Taylor0{Float64} ,nextInputTime::MVector{T,Float64}, x::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T}
  df=0.0
  oldDerDerX=((x[i].coeffs[3])*2.0)
  #@show x
  newDerDerX=(tt.coeffs[2])# 1st der of tt cuz tt itself is derx=f
  #@show newDerDerX
  if elapsed > 0.0
      df=(newDerDerX-oldDerDerX)/(elapsed)
     # println("df=new-old= ",df)
  else
      df= quantum[i]*1e6#*1e12
  end       
  if df!=0.0
      nextInputTime[i]=currentTime+cbrt((abs(quantum[i]/df)))    #df mimics 3rd der 
     # println("usedcbrt")
  else
      if newDerDerX==0 && x[i][1]!=0 #predicted second derivative is 0 & should be not be used to determine nexttime. use 1st der
         #nextInputTime[i]=currentTime+(abs(2*quantum[i] / x[i][1]))  #*2 is not analytic: is just there to increase stepsize
         nextInputTime[i]=currentTime+sqrt(abs(1*quantum[i] / x[i][1]))  #I used the same formulae even with 1st der so that it is fair to other vars
      else
          nextInputTime[i] = Inf
      end
  end
    return nothing
end
#= function computeNextInputTime(::Val{3}, i::Int, currentTime::Float64,elapsed::Float64, tt::Taylor0{Float64} ,nextInputTime::MVector{T,Float64}, x::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T}
  df=0.0
  oldDerDerDerX=((x[i].coeffs[4])*6.0)#@show x
  newDerDerDerX=(tt.coeffs[3])*2.0#@show newDerDerDerX
    if elapsed > 0.0
      df=(newDerDerDerX-oldDerDerDerX)/(elapsed/3) 
    else
      df= quantum[i]*1e6#*1e18   
    end       
   if df!=0.0
    nextInputTime[i]=currentTime+((abs(quantum[i]/df)))^0.25      #df mimics 4th der
   else
      if newDerDerDerX==0 #predicted 3rd derivative is 0 & should be not be used to determine nexttime. use 2nd or 1st der
        if x[i][2]!=0 
            nextInputTime[i]=currentTime+cbrt(abs(6*quantum[i] / x[i][2])) #I used the same formulae even with 2nd der so that it is fair to other vars
          
        elseif x[i][2]==0 && x[i][1]!=0 # second derivative is 0 & should be not be used to determine nexttime. use 1st der
          nextInputTime[i]=currentTime+cbrt(abs(1*quantum[i] / x[i][1]))# #I used the same formulae even with 1st der so that it is fair to other vars
          # println("hllllll")
        else
          nextInputTime[i] = Inf
        end  
      else
        nextInputTime[i] = Inf
      end
   end
    return nothing
end =#

###########################################################################################################################################################""
#= function computeNextEventTime(j::Int,ZCFun::Float64,oldsignValue,currentTime,  nextEventTime, quantum::Vector{Float64})#,printCounter::Vector{Int}) #later specify args
  if oldsignValue[j,1] != sign(ZCFun)
    nextEventTime[j]=currentTime 
  else
    #nextEventTime[j] =currentTime + minPosRoot(ZCFun.coeffs, Val(2)) #Inf  # we can estimate the time. this is important when zc depends only on time   
    nextEventTime[j]=Inf
  end
  oldsignValue[j,1]=sign(ZCFun)#update the values
  oldsignValue[j,2]=ZCFun
end =#

function computeNextEventTime(j::Int,ZCFun::Taylor0{Float64},oldsignValue,currentTime,  nextEventTime, quantum::Vector{Float64})#,printCounter::Vector{Int}) #later specify args
  if oldsignValue[j,1] != sign(ZCFun[0]) || ZCFun[0]==0.0
    nextEventTime[j]=currentTime 
  else
    coef=@SVector [ZCFun[0],ZCFun[1],ZCFun[2]]
    mpr=minPosRoot(coef, Val(2)) 
   #=  x=coef[1]+coef[2]*mpr+coef[3]*mpr*mpr
      @show x =#
    nextEventTime[j] =currentTime + mpr
    #nextEventTime[j] =currentTime + minPosRoot(ZCFun.coeffs, Val(2)) #Inf  # we can estimate the time. this is important when zc depends only on time   
   # nextEventTime[j]=Inf
  # @show currentTime,nextEventTime[j]
    #= if currentTime>=1.132515749331604
      x=coef[1]+coef[2]*mpr+coef[3]*mpr*mpr
      @show x
    end =#
  end
  oldsignValue[j,1]=sign(ZCFun[0])#update the values
  oldsignValue[j,2]=ZCFun[0]
end