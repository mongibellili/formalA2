struct Sol{T}
  order::Int
  savedTimes::Vector{Float64} 
  savedVars::Vector{Array{Taylor0{Float64}}}
  algName::String
  sysName::String
  numSteps ::MVector{T,Int}
  absQ::Float64
  simulStepCount::Int
end

####################################################################################################
function evaluateSol(sol::Sol{T},index::Int,t::Float64)where {T}
  ord=sol.order
  x=Taylor0(zeros(ord + 1), ord) 
  integratorCache=Taylor0(zeros(ord+1),ord)
  for i=2:length(sol[1])#savedTimes after the init time...init time is at index i=1
      if sol[1][i]>t # i-1 is closest lower point
        elapsed=t-sol[1][i-1]  
        x.coeffs .=sol[2][index][i-1].coeffs
       # println("x before= ",x)
        integrateState(Val(ord),x,integratorCache,elapsed)
       # @show t
        #@show sol[1][i]
      #  @show x
        return x#taylor evaluation after small elapsed with the point before (i-1)
      elseif sol[1][i]==t # i-1 is closest lower point
        x.coeffs .=sol[2][index][i].coeffs
      #  @show t
       # @show x
        return x
      end
  end
end
function solInterpolated(sol::Sol{T},step::Float64,ft::Float64)where {T}
  numPoints=length(sol.savedTimes)
  interpTimes=Float64[]
  t=0.0  #later can change to init_time which could be diff than zero
  push!(interpTimes,t)
  while t+step<ft
    t=t+step
    push!(interpTimes,t) 
  end
  push!(interpTimes,ft)
  numInterpPoints=length(interpTimes)
  #display(interpTimes)
  interpValues=Vector{Array{Taylor0{Float64}}}(undef, T)
  for index=1:T
    interpValues[index]=[]
    push!(interpValues[index],sol[2][index][1]) #1st element is the init cond (true value)
  end
  for i=2:numInterpPoints-1
    for index=1:T
      push!(interpValues[index],evaluateSol(sol,index,interpTimes[i]))
    end
  end
  for index=1:T
    push!(interpValues[index],sol[2][index][numPoints]) #last pt @ft
  end
  #(interpTimes,interpValues)
  Sol(sol.order,interpTimes, interpValues,sol.algName,sol.sysName,sol.numSteps,sol.absQ,sol.simulStepCount)
end

function evaluateSimpleSol(sol::Sol,index::Int,t::Float64)
  for i=2:length(sol[1])#savedTimes after the init time...init time is at index i=1
      if sol[1][i]>=t # i-1 is closest lower point
          return sol[2][index][i-1](t-sol[1][i-1])#taylor evaluation after small elapsed with the point before (i-1)
      end
  end
end
function simpleSolInterpolated(sol::Sol,index::Int,step::Float64,ft::Float64)
  numPoints=length(sol.savedTimes)
  interpTimes=[]
  t=0.0  #later can change to init_time which could be diff than zero
  push!(interpTimes,t)
  while t+step<ft
    t=t+step
    push!(interpTimes,t)
    
  end
  push!(interpTimes,ft)
  numInterpPoints=length(interpTimes)
  #display(interpTimes)
  interpValues=[]
  push!(interpValues,sol[2][index][1][0]) #1st element is the init cond (true value)
  #k=2 #k=1 corresponds to initvalue 
  #= for i=2:numPoints#savedTimes after the init time...init time is at index i=1
      if sol[1][i]>interpTimes[k] # i-1 is closest lower point
          push!(interpValues, sol[2][index][i-1](interpTimes[k]-sol[1][i-1]))#taylor evaluation after small elapsed with the point before (i-1)    
          k=k+1  #since no sol[1][i]>ft=interpTimes[last], we are safe to access  interpTimes[last+1]  
      end
  end =#
  for i=2:numInterpPoints-1
    push!(interpValues,evaluateSol(sol,index,interpTimes[i]))
  end
  push!(interpValues,sol[2][index][numPoints][0]) #last pt @ft
  (interpTimes,interpValues)
end
(sol::Sol)(index::Int,t::Float64) = evaluateSol(sol,index,t)
####################################################################################################












####################################################################################################


function getError(sol::Sol,index::Int,f::Function)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
  sumTrueSqr=0.0
  sumDiffSqr=0.0
  relerror=0.0
  if index<=numVars
    for i = 1:numPoints #each point is a taylor
      ft=f(sol.savedTimes[i])
      sumDiffSqr+=(sol.savedVars[index][i].coeffs[1]-ft)*(sol.savedVars[index][i].coeffs[1]-ft)
      sumTrueSqr+=ft*ft
    end
    relerror=sqrt(sumDiffSqr/sumTrueSqr)
  else
    error("the system contains only $numVars variables!")
  end
  return relerror
end
function getIntervalError(sol::Sol,index::Int,f::Function,t0::Float64,tf::Float64)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
  sumTrueSqr=0.0
  sumDiffSqr=0.0
  relerror=0.0
  if index<=numVars
    for i = 1:numPoints #each point is a taylor
      t=sol.savedTimes[i]
      if t0<=t<=tf
        ft=f(t)
        sumDiffSqr+=(sol.savedVars[index][i].coeffs[1]-ft)*(sol.savedVars[index][i].coeffs[1]-ft)
        sumTrueSqr+=ft*ft
      end
    end
    if sumTrueSqr!=0.0
       relerror=sqrt(sumDiffSqr/sumTrueSqr)
    end
  else
    error("the system contains only $numVars variables!")
  end
  return relerror
end
function plotRelativeError(sol::Sol,index::Int,f::Function)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
  if index<=numVars
    temp = []
    tempt = []
    for i = 1:numPoints #each point is a taylor
      ft=f(sol.savedTimes[i])
      if ft>1e-12 || ft < -1e-12
        push!(temp, abs((sol.savedVars[index][i].coeffs[1]-ft)/ft))
        push!(tempt,sol.savedTimes[i])
      end
    end
    display(plot(tempt, temp,title="RelError:(S-T)/T for x$(index)_$(sol.absQ)",label="$(sol.algName)"#= ,xlims=(80,200),ylims=(0.0,0.0002) =#) )
    println("press enter to exit")
    readline() 
  else
    error("the system contains only $numVars variables!")
  end
end
function stackplotRelativeError(sol::Sol,index::Int,f::Function)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
  if index<=numVars
    temp = []
    tempt = []
    for i = 1:numPoints #each point is a taylor
      ft=f(sol.savedTimes[i])
      if ft>1e-12 || ft < -1e-12
        push!(temp, abs((sol.savedVars[index][i].coeffs[1]-ft)/ft))
        push!(tempt,sol.savedTimes[i])
      end
    end
    display(plot!(tempt, temp,title="RelError:(S-T)/T for x$(index)_$(sol.absQ)",label="$(sol.algName)"#= ,xlims=(80,200),ylims=(0.0,0.0002) =#) )
    println("press enter to exit")
    readline() 
  else
    error("the system contains only $numVars variables!")
  end
end
function plot_save_RelativeError(sol::Sol,index::Int,f::Function)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)

  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
  if index<=numVars
    temp = []
    tempt = []
    for i = 1:numPoints #each point is a taylor
      ft=f(sol.savedTimes[i])
      if ft>1e-12 || ft < -1e-12
        push!(temp, abs((sol.savedVars[index][i].coeffs[1]-ft)/ft))
        push!(tempt,sol.savedTimes[i])
      end
    end
  # display(plot!(tempt, temp,title="RelError:(S-T)/T for x$index",label="$(sol.algName)")) 
  p1=plot(tempt, temp,title="RelError:(S-T)/T for x$(index)_$(sol.absQ)",label="$(sol.algName)"#= ,xlims=(80,200),ylims=(0.0,0.0002) =#) 
  savefig(p1, "relError_$(sol.sysName)_$(sol.algName)_$(sol.absQ)_x$(index)_$(timestamp).png")
  display(p1)
 println("press enter to exit")
 readline() 
  else
    error("the system contains only $numVars variables!")
  end
end
function stackplot_save_RelativeError(sol::Sol,index::Int,f::Function)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)

  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
  if index<=numVars
    temp = []
    tempt = []
    for i = 1:numPoints #each point is a taylor
      ft=f(sol.savedTimes[i])
      if ft>1e-12 || ft < -1e-12
        push!(temp, abs((sol.savedVars[index][i].coeffs[1]-ft)/ft))
        push!(tempt,sol.savedTimes[i])
      end
    end
  # display(plot!(tempt, temp,title="RelError:(S-T)/T for x$index",label="$(sol.algName)")) 
  p1=plot(tempt, temp,title="RelError:(S-T)/T for  x$(index)_$(sol.absQ)",label="$(sol.algName)"#= ,xlims=(80,200),ylims=(0.0,0.0002) =#) 
  savefig(p1, "relError_$(sol.sysName)_$(sol.algName)_$(sol.absQ)_x$(index)_$(timestamp).png")
  display(p1)
 println("press enter to exit")
 readline() 
  else
    error("the system contains only $numVars variables!")
  end
end
function saveRelativeError(sol::Sol,index::Int,f::Function)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)

  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
  if index<=numVars
    temp = []
    tempt = []
    for i = 1:numPoints #each point is a taylor
      ft=f(sol.savedTimes[i])
      if ft>1e-12 || ft < -1e-12
        push!(temp, abs((sol.savedVars[index][i].coeffs[1]-ft)/ft))
        push!(tempt,sol.savedTimes[i])
      end
    end
  # display(plot!(tempt, temp,title="RelError:(S-T)/T for x$index",label="$(sol.algName)")) 
  p1=plot(tempt, temp,title="RelError:(S-T)/T for x$(index)_$(sol.absQ)",label="$(sol.algName)"#= ,xlims=(80,200),ylims=(0.0,0.0002) =#) 
  savefig(p1, "relError_$(sol.sysName)_$(sol.algName)_$(sol.absQ)_x$(index)_$(timestamp).png")
  else
    error("the system contains only $numVars variables!")
  end
end
function stacksaveRelativeError(sol::Sol,index::Int,f::Function)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)

  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
  if index<=numVars
    temp = []
    tempt = []
    for i = 1:numPoints #each point is a taylor
      ft=f(sol.savedTimes[i])
      if ft>1e-12 || ft < -1e-12
        push!(temp, abs((sol.savedVars[index][i].coeffs[1]-ft)/ft))
        push!(tempt,sol.savedTimes[i])
      end
    end
  # display(plot!(tempt, temp,title="RelError:(S-T)/T for x$index",label="$(sol.algName)")) 
  p1=plot!(tempt, temp,title="RelError:(S-T)/T for x$(index)_$(sol.absQ)",label="$(sol.algName)"#= ,xlims=(80,200),ylims=(0.0,0.0002) =#) 
  savefig(p1, "relError_$(sol.sysName)_$(sol.algName)_$(sol.absQ)_x$(index)_$(timestamp).png")
  else
    error("the system contains only $numVars variables!")
  end
end
function plotAbsoluteError(sol::Sol,index::Int,f::Function)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
  if index<=numVars
    temp = []
   # tempt = []
    for i = 1:numPoints #each point is a taylor
      ft=f(sol.savedTimes[i])
      #if ft>1e-2 || ft < -1e-2
        push!(temp, abs((sol.savedVars[index][i].coeffs[1]-ft)))
       # push!(tempt,sol.savedTimes[i])
     # end
    end
   #display(plot!(sol.savedTimes, temp,title="AbsError:(S-T) for x$index",label="$(sol.algName)")) 
  display(plot(sol.savedTimes, temp,title="AbsError:(S-T) for x$(index)_$(sol.absQ)",label="$(sol.algName)"#= ,xlims=(80,200),ylims=(0.0,0.0002) =#))
  
    println("press enter to exit")
    readline() 
  else
    error("the system contains only $numVars variables!")
  end
end

function stackplotAbsoluteError(sol::Sol,index::Int,f::Function)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
  if index<=numVars
    temp = []
   # tempt = []
    for i = 1:numPoints #each point is a taylor
      ft=f(sol.savedTimes[i])
      #if ft>1e-2 || ft < -1e-2
        push!(temp, abs((sol.savedVars[index][i].coeffs[1]-ft)))
       # push!(tempt,sol.savedTimes[i])
     # end
    end
   #display(plot!(sol.savedTimes, temp,title="AbsError:(S-T) for x$index",label="$(sol.algName)")) 
   display(plot!(sol.savedTimes, temp,title="AbsError:(S-T) for x$(index)_$(sol.absQ)",label="$(sol.algName)"#= ,xlims=(80,200),ylims=(0.0,0.0002) =#))
  
    println("press enter to exit")
    readline() 
  else
    error("the system contains only $numVars variables!")
  end
end

function plotsaveAbsoluteError(sol::Sol,index::Int,f::Function)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)

  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
  if index<=numVars
    temp = []
   # tempt = []
    for i = 1:numPoints #each point is a taylor
      ft=f(sol.savedTimes[i])
      #if ft>1e-2 || ft < -1e-2
        push!(temp, abs((sol.savedVars[index][i].coeffs[1]-ft)))
       # push!(tempt,sol.savedTimes[i])
     # end
    end
   #display(plot!(sol.savedTimes, temp,title="AbsError:(S-T) for x$index",label="$(sol.algName)")) 
   p1=plot(sol.savedTimes, temp,title="AbsError:(S-T) for x$(index)_$(sol.absQ)",label="$(sol.algName)"#= ,xlims=(80,200),ylims=(0.0,0.0002) =#)
   savefig(p1, "absError_$(sol.sysName)_$(sol.algName)_$(sol.absQ)_x$(index)_$(timestamp).png")
     display(p1)
    println("press enter to exit")
    readline() 
  else
    error("the system contains only $numVars variables!")
  end
end
function stackplotsaveAbsoluteError(sol::Sol,index::Int,f::Function)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
  p1=nothing
  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
  if index<=numVars
    temp = []
   # tempt = []
    for i = 1:numPoints #each point is a taylor
      ft=f(sol.savedTimes[i])
      #if ft>1e-2 || ft < -1e-2
        push!(temp, abs((sol.savedVars[index][i].coeffs[1]-ft)))
       # push!(tempt,sol.savedTimes[i])
     # end
    end
   #display(plot!(sol.savedTimes, temp,title="AbsError:(S-T) for x$index",label="$(sol.algName)")) 
   p1=plot!(sol.savedTimes, temp,title="AbsError:(S-T) for x$(index)_$(sol.absQ)",label="$(sol.algName)"#= ,xlims=(80,200),ylims=(0.0,0.0002) =#)
   savefig(p1, "absError_$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(timestamp).png")
     display(p1)
    println("press enter to exit")
    readline() 
  else
    error("the system contains only $numVars variables!")
  end
end
function saveAbsoluteError(sol::Sol,index::Int,f::Function)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
 
  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
  if index<=numVars
    temp = []
   # tempt = []
    for i = 1:numPoints #each point is a taylor
      ft=f(sol.savedTimes[i])
      #if ft>1e-2 || ft < -1e-2
        push!(temp, abs((sol.savedVars[index][i].coeffs[1]-ft)))
       # push!(tempt,sol.savedTimes[i])
     # end
    end
   #display(plot!(sol.savedTimes, temp,title="AbsError:(S-T) for x$index",label="$(sol.algName)")) 
   p1=plot(sol.savedTimes, temp,title="AbsError:(S-T) for x$(index)_$(sol.absQ)",label="$(sol.algName)"#= ,xlims=(80,200),ylims=(0.0,0.0002) =#)
   savefig(p1, "absError_$(sol.sysName)_$(sol.algName)_$(sol.absQ)_x$(index)_$(timestamp).png")
     
  else
    error("the system contains only $numVars variables!")
  end
end
function saveZoomedAbsoluteError(sol::Sol,index::Int,f::Function)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
 
  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
  if index<=numVars
    temp = []
   # tempt = []
    for i = 1:numPoints #each point is a taylor
      ft=f(sol.savedTimes[i])
      #if ft>1e-2 || ft < -1e-2
        push!(temp, abs((sol.savedVars[index][i].coeffs[1]-ft)))
       # push!(tempt,sol.savedTimes[i])
     # end
    end
   #display(plot!(sol.savedTimes, temp,title="AbsError:(S-T) for x$index",label="$(sol.algName)")) 
   p1=plot(sol.savedTimes, temp,title="AbsError:(S-T) for x$(index)_$(sol.absQ)",label="$(sol.algName)"#= ,xlims=(80,200)=#,ylims=(0.0,5e-6) )
   savefig(p1, "absError_$(sol.sysName)_$(sol.algName)_$(sol.absQ)_x$(index)_$(timestamp).png")
     
  else
    error("the system contains only $numVars variables!")
  end
end
function stacksaveAbsoluteError(sol::Sol,index::Int,f::Function)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)

  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
  if index<=numVars
    temp = []
   # tempt = []
    for i = 1:numPoints #each point is a taylor
      ft=f(sol.savedTimes[i])
      #if ft>1e-2 || ft < -1e-2
        push!(temp, abs((sol.savedVars[index][i].coeffs[1]-ft)))
       # push!(tempt,sol.savedTimes[i])
     # end
    end
   #display(plot!(sol.savedTimes, temp,title="AbsError:(S-T) for x$index",label="$(sol.algName)")) 
   p1=plot!(sol.savedTimes, temp,title="AbsError:(S-T) for x$(index)_$(sol.absQ)",label="$(sol.algName)"#= ,xlims=(80,200),ylims=(0.0,0.0002) =#)
   savefig(p1, "absError_$(sol.sysName)_$(sol.algName)_$(sol.absQ)_x$(index)_$(timestamp).png")
    
  else
    error("the system contains only $numVars variables!")
  end
end
function plotCumulativeSquaredRelativeError(sol::Sol,index::Int,f::Function)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
  sumTrueSqr=0.0
  sumDiffSqr=0.0
  
  if index<=numVars
    temp = []
    for i = 1:numPoints #each point is a taylor
      ft=f(sol.savedTimes[i])
      sumDiffSqr+=(sol.savedVars[index][i].coeffs[1]-ft)*(sol.savedVars[index][i].coeffs[1]-ft)
      sumTrueSqr+=ft*ft
        push!(temp, sqrt(sumDiffSqr/sumTrueSqr))
    end
   display(plot!(sol.savedTimes, temp,title="Error:sqrt(∑(S-T)^2/∑T^2) for x$index",label="$(sol.algName)")) 
  else
    error("the system contains only $numVars variables!")
  end
  println("press enter to exit")
  readline()
end
function plotMSE(sol::Sol,index::Int,f::Function)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
 # sumTrueSqr=0.0
  sumDiffSqr=0.0
  
  if index<=numVars
    temp = []
    for i = 1:numPoints #each point is a taylor
      ft=f(sol.savedTimes[i])
      sumDiffSqr+=(sol.savedVars[index][i].coeffs[1]-ft)*(sol.savedVars[index][i].coeffs[1]-ft)
     # sumTrueSqr+=ft*ft
        push!(temp, (sumDiffSqr/i))
    end
   display(plot!(sol.savedTimes, temp,title="Error:(∑(S-T)^2/i) for x$index",label="$(sol.algName)")) 
  else
    error("the system contains only $numVars variables!")
  end
  println("press enter to exit")
  readline()
end
####################################################################################################










function plotSol(sol::Sol)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
  p1=plot()
  for k=1:numVars
    temp = []
    for i = 1:numPoints #each point is a taylor
        push!(temp, sol.savedVars[k][i].coeffs[1])
    end
    display(plot!(p1,sol.savedTimes, temp,marker=(:circle),markersize=3,title="$(sol.algName)  _$(sol.absQ)_$(sol.simulStepCount)",label="x$k:$(sol.numSteps[k]) steps"#= ,xlims=(877.5,877.8),ylims=(-0.01,0.03) =#))
  end
    println("press enter to exit")
    readline() 
end
function stackplotSol(sol::Sol)
    numPoints=length(sol.savedTimes)
    numVars=length(sol.savedVars)
   # p1=nothing;p2=nothing
    for k=1:numVars
      temp = []
      for i = 1:numPoints #each point is a taylor
          push!(temp, sol.savedVars[k][i].coeffs[1])
      end
     
      display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.sysName)_$(sol.absQ)",label="$(sol.algName)_x$k:$(sol.numSteps[k]) steps"#= ,xlims=(877.5,877.8),ylims=(-0.01,0.03) =#))
     
      
   #  display(plot!(sol.savedTimes,marker=(:circle), temp,title="$(sol.algName)",label="$(sol.algName):-x$k"#= ,xlims=(435,530),ylims=(0,1) =#))
   #display(plot!(sol.savedTimes,marker=(:circle), temp,title="$(sol.algName)",label="$(sol.algName):-x$k",xlims=(62.2,62.5),ylims=(9.8,10.35)))  #both medium zoom
  # display(plot!(sol.savedTimes,marker=(:circle), temp,title="$(sol.algName)",label="$(sol.algName):-x$k",xlims=(63.2,63.4),ylims=(9.9,10.3))) #mliqss2 zoom
   #display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(435,800),ylims=(-0.02,0.1)))
   #display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(800,2000),ylims=(-0.02,0.1)))#+1e-19
  # display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(877.5,877.8),ylims=(-0.01,0.03)))#++1e-19
 #=   println("press enter to exit")
   readline() =#
   #x1
 #  display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(918,920),ylims=(20.1889,20.1890))) #+eps=1e-6
 #  display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(700,2000),ylims=(20.17,20.21)))#++1e-19 
 #  display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(842,850),ylims=(20.1873,20.1874)))#++1e-19
 #display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(847.35,847.65),ylims=(20.1873,20.1874)))#+mliqss2+1e-19
 #display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(877.5,877.8),ylims=(20.19475,20.19480)))#+liqss2+1e-19
   #  display(plot!(sol.savedTimes, temp,title="System1-qss2",label="x$k",xlims=(0,0.6),ylims=(-0.4,1))) 
    #=  if k%2==1#k=1 sys5 x1
      display(plot!(sol.savedTimes,temp#= ,line=(:dot,3) =#,marker=(:circle)#= marker=([:cross :d]) =#,title="$(sol.algName)",label="$(sol.algName): x1",xlims=(80,200),ylims=(0.7495,0.7505))) # system5 against true solution
     end =#
    #=  if k%2==0#k=2 sys5 x2
      display(plot!(sol.savedTimes,temp#= ,line=(:dot,3) =#,marker=(:circle)#= marker=([:cross :d]) =#,title="$(sol.algName)",label="$(sol.algName): x2",xlims=(80,200),ylims=(3.999,4.0020))) # system5 against true solution
     end =#
    #=  if k%2==0#k=2 sys3 x2
      display(plot!(sol.savedTimes,temp,#= seriestype = :scatter =##= line=(:dot,3), =#marker=(:circle),title="$(sol.algName)",label="$(sol.algName): x2",xlims=(7,50),ylims=(0.6998,0.7003)))
     # display(plot!(sol.savedTimes,temp,#= seriestype = :scatter =#line=(:dot,3),title="$(sol.algName)",label="$(sol.algName): x2",xlims=(2,20),ylims=(0.6991,0.701)))
     # display(plot!(sol.savedTimes,temp,#= seriestype = :scatter =#line=(:dot,3),title="$(sol.algName)",label="$(sol.algName): x2",xlims=(2,15),ylims=(0.698,0.81)))
     end =#
    

    # display(plot!(sol.savedTimes, temp,seriestype = :scatter,label="x$k",xlims=(10,12),ylims=(-0.5,0.5)))  
    #------------ twoVarSys1----------------
    #-----medium zoom---
    #display(plot!(sol.savedTimes, temp,label="x$k",xlims=(20,50),ylims=(-0.00003,0.00003))) 
      #------------ twoVarSys13----------------
      #mliqss1 adequate
      #= if k%2==0
      display(plot!(sol.savedTimes,label="x$k",xlims=(10,50),ylims=(-0.02,0.1))) 
      else
        display(plot!(sol.savedTimes,label="x$k",line=(:dot,4),xlims=(10,50),ylims=(-0.02,0.1))) 
      end =#
      #liqss2 adequate
      #= if k%2==0
        display(plot!(sol.savedTimes,label="x$k",xlims=(10,50),ylims=(-0.02,1.1))) #1.1 for small delta, 3.1 for a larger delta
        else
          display(plot!(sol.savedTimes,label="x$k",line=(:dot,4),xlims=(10,50),ylims=(-0.02,1.1))) 
        end =#
      # display(plot!(sol.savedTimes, temp,label="x$k",xlims=(10.0,10.0000001),ylims=(0.0318803,0.0318810))) 
            
          
          # display(plot!(sol.savedTimes, temp,label="x$k",xlims=(10.0,10.0000001),ylims=(0.0318803,0.0318810))) 
        # display(plot!(sol.savedTimes, temp,label="x$k",xlims=(8.8,11.2),ylims=(-0.5,0.5))) 
          # display(plot!(sol.savedTimes, temp,label="x$k",xlims=(10.0,10.0000001),ylims=(0.0318803,0.0318810))) 
      
    end
   
    
      # display(plot(p1, p2, layout = 2))
      println("press enter to exit")
      readline() 
     
end

function plot_save_Sol(sol::Sol)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
  p1=plot()#;p2=nothing
  #timestamp=string(round(now(),Minute(1)))# remove colon file name complains
  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
  for k=1:numVars
    temp = []
    for i = 1:numPoints #each point is a taylor
        push!(temp, sol.savedVars[k][i].coeffs[1])
    end
    p1=plot!(p1,sol.savedTimes, temp,marker=(:circle),title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)",label="x$k:$(sol.numSteps[k]) steps"#= ,xlims=(877.5,877.8),ylims=(-0.01,0.03) =#)
  end
  savefig(p1, "plot_$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(timestamp).png")
  display(p1)
  println("press enter to exit")
  readline()  
end
function stackplot_save_Sol(sol::Sol)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
  p1=nothing#;p2=nothing
  #timestamp=string(round(now(),Minute(1)))# remove colon file name complains
  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
  for k=1:numVars
    temp = []
    for i = 1:numPoints #each point is a taylor
        push!(temp, sol.savedVars[k][i].coeffs[1])
    end
    p1=plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.sysName)_$(sol.absQ)",label="_$(sol.algName)_x$(k):$(sol.numSteps[k]) steps"#= ,xlims=(877.5,877.8),ylims=(-0.01,0.03) =#)
  end
  savefig(p1, "plot_$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(timestamp).png")
  display(p1)
  println("press enter to exit")
  readline()    
end

function save_Sol(sol::Sol)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
  p1=plot()#;p2=nothing
  #timestamp=string(round(now(),Minute(1)))# remove colon file name complains
  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
  for k=1:numVars
    temp = []
    for i = 1:numPoints #each point is a taylor
        push!(temp, sol.savedVars[k][i].coeffs[1])
    end
    p1=plot!(p1,sol.savedTimes, temp,marker=(:circle),title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.simulStepCount)",label="x$k:$(sol.numSteps[k]) steps"#= ,xlims=(877.5,877.8),ylims=(-0.01,0.03) =#)
    #p1=plot!(p1,sol.savedTimes, temp,title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.simulStepCount)",label="x$k:$(sol.numSteps[k]) steps"#= ,xlims=(877.5,877.8),ylims=(-0.01,0.03) =#)
  end
  savefig(p1, "plot_$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(timestamp).png")
end
function stack_save_Sol(sol::Sol)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
  p1=nothing#;p2=nothing
  #timestamp=string(round(now(),Minute(1)))# remove colon file name complains
  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
  for k=1:numVars
    temp = []
    for i = 1:numPoints #each point is a taylor
        push!(temp, sol.savedVars[k][i].coeffs[1])
    end
    p1=plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.sysName)_$(sol.absQ)",label="_$(sol.algName)_x$k:$(sol.numSteps[k]) steps"#= ,xlims=(877.5,877.8),ylims=(-0.01,0.03) =#)
  end
  savefig(p1, "plot_$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(timestamp).png")  
end




function plot_save_SolVars(sol::Sol)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
  p1=nothing;p2=nothing
  #timestamp=string(round(now(),Minute(1)))# remove colon file name complains
  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
  for k=1:numVars
    temp = []
    for i = 1:numPoints #each point is a taylor
        push!(temp, sol.savedVars[k][i].coeffs[1])
    end
    if k%2==0#k=2
    p1=plot(sol.savedTimes, temp,marker=(:circle),title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)",label="x$k:$(sol.numSteps[k]) steps"#= ,xlims=(877.5,877.8),ylims=(-0.01,0.03) =#)
    
    savefig(p1, "plot_$(sol.sysName)_x$(k)_$(sol.algName)_$(sol.absQ)_$(timestamp).png")
    #savefig(p1, "liqx$k.png")
    else
    p2=plot(sol.savedTimes, temp,marker=(:circle),title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)",label="x$k:$(sol.numSteps[k]) steps"#= ,xlims=(877.5,877.8),ylims=(20.19475,20.19480) =#)
    savefig(p2, "plot_$(sol.sysName)_x$(k)_$(sol.algName)_$(sol.absQ)_$(timestamp).png")
    end
    
      #  display(plot!(sol.savedTimes,marker=(:circle), temp,title="$(sol.algName)",label="$(sol.algName):-x$k"#= ,xlims=(435,530),ylims=(0,1) =#))
      #display(plot!(sol.savedTimes,marker=(:circle), temp,title="$(sol.algName)",label="$(sol.algName):-x$k",xlims=(62.2,62.5),ylims=(9.8,10.35)))  #both medium zoom
      # display(plot!(sol.savedTimes,marker=(:circle), temp,title="$(sol.algName)",label="$(sol.algName):-x$k",xlims=(63.2,63.4),ylims=(9.9,10.3))) #mliqss2 zoom
      #display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(435,800),ylims=(-0.02,0.1)))
      #display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(800,2000),ylims=(-0.02,0.1)))#+1e-19
      # display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(877.5,877.8),ylims=(-0.01,0.03)))#++1e-19
      #=   println("press enter to exit")
      readline() =#
      #x1
      #  display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(918,920),ylims=(20.1889,20.1890))) #+eps=1e-6
      #  display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(700,2000),ylims=(20.17,20.21)))#++1e-19 
      #  display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(842,850),ylims=(20.1873,20.1874)))#++1e-19
      #display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(847.35,847.65),ylims=(20.1873,20.1874)))#+mliqss2+1e-19
      #display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(877.5,877.8),ylims=(20.19475,20.19480)))#+liqss2+1e-19
      #  display(plot!(sol.savedTimes, temp,title="System1-qss2",label="x$k",xlims=(0,0.6),ylims=(-0.4,1))) 
        #=  if k%2==1#k=1 sys5 x1
          display(plot!(sol.savedTimes,temp#= ,line=(:dot,3) =#,marker=(:circle)#= marker=([:cross :d]) =#,title="$(sol.algName)",label="$(sol.algName): x1",xlims=(80,200),ylims=(0.7495,0.7505))) # system5 against true solution
        end =#
        #=  if k%2==0#k=2 sys5 x2
          display(plot!(sol.savedTimes,temp#= ,line=(:dot,3) =#,marker=(:circle)#= marker=([:cross :d]) =#,title="$(sol.algName)",label="$(sol.algName): x2",xlims=(80,200),ylims=(3.999,4.0020))) # system5 against true solution
        end =#
        #=  if k%2==0#k=2 sys3 x2
          display(plot!(sol.savedTimes,temp,#= seriestype = :scatter =##= line=(:dot,3), =#marker=(:circle),title="$(sol.algName)",label="$(sol.algName): x2",xlims=(7,50),ylims=(0.6998,0.7003)))
        # display(plot!(sol.savedTimes,temp,#= seriestype = :scatter =#line=(:dot,3),title="$(sol.algName)",label="$(sol.algName): x2",xlims=(2,20),ylims=(0.6991,0.701)))
        # display(plot!(sol.savedTimes,temp,#= seriestype = :scatter =#line=(:dot,3),title="$(sol.algName)",label="$(sol.algName): x2",xlims=(2,15),ylims=(0.698,0.81)))
        end =#
        

        # display(plot!(sol.savedTimes, temp,seriestype = :scatter,label="x$k",xlims=(10,12),ylims=(-0.5,0.5)))  
        #------------ twoVarSys1----------------
        #-----medium zoom---
        #display(plot!(sol.savedTimes, temp,label="x$k",xlims=(20,50),ylims=(-0.00003,0.00003))) 
          #------------ twoVarSys13----------------
          #mliqss1 adequate
          #= if k%2==0
          display(plot!(sol.savedTimes,label="x$k",xlims=(10,50),ylims=(-0.02,0.1))) 
          else
            display(plot!(sol.savedTimes,label="x$k",line=(:dot,4),xlims=(10,50),ylims=(-0.02,0.1))) 
          end =#
          #liqss2 adequate
          #= if k%2==0
            display(plot!(sol.savedTimes,label="x$k",xlims=(10,50),ylims=(-0.02,1.1))) #1.1 for small delta, 3.1 for a larger delta
            else
              display(plot!(sol.savedTimes,label="x$k",line=(:dot,4),xlims=(10,50),ylims=(-0.02,1.1))) 
            end =#
          # display(plot!(sol.savedTimes, temp,label="x$k",xlims=(10.0,10.0000001),ylims=(0.0318803,0.0318810))) 
                
              
              # display(plot!(sol.savedTimes, temp,label="x$k",xlims=(10.0,10.0000001),ylims=(0.0318803,0.0318810))) 
            # display(plot!(sol.savedTimes, temp,label="x$k",xlims=(8.8,11.2),ylims=(-0.5,0.5))) 
              # display(plot!(sol.savedTimes, temp,label="x$k",xlims=(10.0,10.0000001),ylims=(0.0318803,0.0318810))) 
    
  end
 
  
    #=  display(plot(p1, p2, layout = 2))
    println("press enter to exit")
    readline()  =#
   
end
function plotSol_Der1(sol::Sol)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
  for k=1:numVars
    temp = []
    for i = 1:numPoints #each point is a taylor
        push!(temp, sol.savedVars[k][i].coeffs[1])
    end
   display(plot!(sol.savedTimes, temp,label="x$k")) 
   temp2 = []
    for i = 1:numPoints #each point is a taylor
        push!(temp2, sol.savedVars[k][i].coeffs[2])
    end
   display(plot!(sol.savedTimes, temp2,title="$(sol.algName)",label="dx$k")) 
    
  end
    println("press enter to exit")
    readline()
end
function plotSol(savedTimes::Vector{Float64} ,  savedVars::Vector{Array{Taylor0{Float64}}})
    numPoints=length(savedTimes)
    numVars=length(savedVars)
    for k=1:numVars
      temp = []
      for i = 1:numPoints #each point is a taylor
          push!(temp, savedVars[k][i].coeffs[1])
      end
      display(plot!(savedTimes, temp,label="x$k"))
    end
      println("press enter to exit")
      readline()
end
function plotSol(savedTimes::Vector{Float64} , savedVar::Vector{Taylor0{Float64}})
      numPoints=length(savedTimes)
      temp = []
      for i = 1:numPoints
          push!(temp, savedVar[i].coeffs[1])
      end
      display(plot!(savedTimes, temp))
      println("press enter to exit")
      readline()
end


  
function getindex(s::Sol, i::Int64)
      if i==1
         return s.savedTimes
      elseif i==2
         return s.savedVars
      else
         error("sol has 2 attributes: time and states")
      end
end

function save_SolVar(sol::Sol,k::Int)
  numPoints=length(sol.savedTimes)
  
  p1=plot()#;p2=nothing
  #timestamp=string(round(now(),Minute(1)))# remove colon file name complains
  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
  
    temp = []
    for i = 1:numPoints #each point is a taylor
        push!(temp, sol.savedVars[k][i].coeffs[1])
    end
    p1=plot!(p1,sol.savedTimes, temp,marker=(:circle),title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.simulStepCount)",label="x$k:$(sol.numSteps[k]) steps"#= ,xlims=(877.5,877.8),ylims=(-0.01,0.03) =#)
    #p1=plot!(p1,sol.savedTimes, temp,title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.simulStepCount)",label="x$k:$(sol.numSteps[k]) steps"#= ,xlims=(877.5,877.8),ylims=(-0.01,0.03) =#)
  
  savefig(p1, "plot_$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(timestamp).png")

end