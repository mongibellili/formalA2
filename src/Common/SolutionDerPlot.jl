

#if you want to put a note and lims in the graph
function plot_SolDer(sol::Sol{T,O},xvars::Int...;note=" "::String,xlims=(0.0,0.0)::Tuple{Float64, Float64},ylims=(0.0,0.0)::Tuple{Float64, Float64},legend=:true::Bool) where{T,O}
  p1=plot()
 
  #= if sol.algName=="nmliqss1"
    mycolor=:green
    stle=:dash
  else
    mycolor=:blue
  end =#
  if xvars!=()
    for k in xvars
      if k==1
        stle=:dash
        sze=2
      elseif k==2
        stle=:solid
        sze=4
      elseif k==3
        stle=:dot
        sze=3
      elseif k==4
        stle=:dashdot
        sze=2
      elseif k==5
        stle=:dashdotdot
        sze=2
      else
        sze=1
        stle=:solid
      end
     # p1=plot!(p1,sol.savedTimes[k], sol.savedVarsQ[k],line=(1,mycolor,:dash),marker=(:star),label="q$k $(sol.algName)")
     # p1=plot!(p1,sol.savedTimes[k], sol.savedDers[k]#= ,marker=(:circle) =#,markersize=2,label="x$k ",legend=:bottomright)
      p1=plot!(p1,sol.savedTimes[k], sol.savedDers[k],line=(sze,stle),marker=(:circle),label="dx$k $(sol.numSteps[k])"#= ,legend=:right =#)
    end
  else
    for k=1:T
      if k==1
        mycolor=:red
      else
        mycolor=:purple
      end
     # p1=plot!(p1,sol.savedTimes[k], sol.savedVarsQ[k],line=(1,mycolor,:dash),marker=(:star),label="q$k $(sol.algName)")
      p1=plot!(p1,sol.savedTimes[k], sol.savedDers[k],marker=(:circle),markersize=2,label="dx$k $(sol.numSteps[k])"#= ,legend=:false =#)
    end
  end
  if xlims!=(0.0,0.0) && ylims!=(0.0,0.0) 
    p1=plot!(p1,title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.totalSteps)_$(sol.simulStepCount)_$(sol.evCount) \n $note", xlims=xlims ,ylims=ylims)
  elseif xlims!=(0.0,0.0) && ylims==(0.0,0.0) 
    p1=plot!(p1,title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.totalSteps)_$(sol.simulStepCount)_$(sol.evCount)  \n $note", xlims=xlims #= ,ylims=ylims =#)
  elseif xlims==(0.0,0.0) && ylims!=(0.0,0.0) 
    p1=plot!(p1,title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.totalSteps)_$(sol.simulStepCount)_$(sol.evCount)  \n $note"#= , xlims=xlims  =#,ylims=ylims)
  else
    p1=plot!(p1, title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.totalSteps)_$(sol.simulStepCount)_$(sol.evCount)  \n $note",legend=legend)
  end
  p1
  #savefig(p1, "plot_$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(note)_ft_$(sol.ft)_$(timestamp).png")
end

function save_SolDer(sol::Sol{T,O},xvars::Int...;note=" "::String,xlims=(0.0,0.0)::Tuple{Float64, Float64},ylims=(0.0,0.0)::Tuple{Float64, Float64},legend=:true::Bool) where{T,O}
  
  p1= plot_SolDer(sol,xvars...;note=note,xlims=xlims,ylims=ylims,legend=legend)
  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
  savefig(p1, "plot_$(sol.sysName)_$(sol.algName)_$(xvars)_$(sol.absQ)_$(note)_ftt_$(sol.ft)_$(timestamp).png")
end

#= #for debug to be deleted later: simultaneous steps plot
function save_SimulSol(sol::Sol{T,O},xvars::Int...;note=" "::String,xlims=(0.0,0.0)::Tuple{Float64, Float64},ylims=(0.0,0.0)::Tuple{Float64, Float64}) where{T,O}
  p1=plot()
  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
 
  if xvars[1]!=()
    for k in xvars
      p1=plot!(p1,sol.savedTimes[k], sol.savedDers[k],marker=(:circle),markersize=2,label="x$k $(sol.numSteps[k])")
      p1=plot!(p1,sol.simulStepsTimes[k], sol.simulStepsVals[k],marker=(:circle),markersize=4,label="x$k ")
    end
  else
    for k=1:T
      p1=plot!(p1,sol.savedTimes[k], sol.savedDers[k],marker=(:circle),markersize=2,label="x$k $(sol.numSteps[k])")
      p1=plot!(p1,sol.simulStepsTimes[k], sol.simulStepsVals[k],marker=(:circle),markersize=4,label="x$k ")
    end
  end
  if xlims!=(0.0,0.0) && ylims!=(0.0,0.0) 
    p1=plot!(p1,title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.totalSteps)_$(sol.simulStepCount) \n $note", xlims=xlims ,ylims=ylims)
  elseif xlims!=(0.0,0.0) && ylims==(0.0,0.0) 
    p1=plot!(p1,title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.totalSteps)_$(sol.simulStepCount) \n $note", xlims=xlims #= ,ylims=ylims =#)
  elseif xlims==(0.0,0.0) && ylims!=(0.0,0.0) 
    p1=plot!(p1,title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.totalSteps)_$(sol.simulStepCount) \n $note"#= , xlims=xlims  =#,ylims=ylims)
  else
    p1=plot!(p1, title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.totalSteps)_$(sol.simulStepCount) \n $note")
  end
  savefig(p1, "plot_$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(note)_ft_$(sol.ft)_$(timestamp).png")
end




function getPlot(sol::Sol{T,O}) where{T,O}
  p1=plot(title="$(sol.sysName)")#;p2=nothing
  for k=1:T
    p1=plot!(p1,sol.savedTimes[k], sol.savedDers[k],marker=(:circle),markersize=2,label="x$k $(sol.numSteps[k]) ")
  end
  p1
end
function getPlot(sol::Sol{T,O},k::Int) where{T,O}
  p1=plot!(sol.savedTimes[k], sol.savedDers[k],marker=(:circle),markersize=2,#= title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.simulStepCount) ", =#label="x$index ")
end
function getPlot!(sol::Sol{T,O}) where{T,O}
  p1=plot!(title="$(sol.sysName)")
  if sol.algName=="nmliqss1"
    mycolor=:red
  else
    mycolor=:blue
  end
  for k=1:T
   
    #p1=plot!(p1,sol.savedTimes[k], sol.savedVarsQ[k],line=(1,mycolor,:dash),label="q$k $(sol.algName)")
    p1=plot!(p1,sol.savedTimes[k], sol.savedDers[k],line=(1,mycolor),marker=(:circle),markersize=2,label="x$k $(sol.numSteps[k])$(sol.algName)_$(sol.totalSteps)_$(sol.stepsaftersimul)_$(sol.simulStepCount)")
  end
  p1
end
function getPlot!(sol::Sol{T,O},k::Int) where{T,O} 
  p1=plot!(title="$(sol.sysName)")
    if sol.algName=="nmliqss1"
      mycolor=:red;linsty=:dash
    else
      mycolor=:blue;linsty=:dot
    end
   # p1=plot!(p1,sol.savedTimes[k], sol.savedVarsQ[k],line=(1,mycolor,:dash),marker=(:star),label="q$k $(sol.algName)")
    p1=plot!(p1,sol.savedTimes[k], sol.savedDers[k],line=(1,mycolor),marker=(:circle),markersize=2,label="x$k $(sol.numSteps[k])$(sol.algName)_$(sol.totalSteps)_$(sol.simulStepCount)")
end
 =#

#= function plotSol(sol::Sol{T,O}) where{T,O}
  numPoints=length(sol.savedTimes)
  #numVars=length(sol.savedDers)
  p1=plot()
  for k=1:T
    temp = []
    for i = 1:numPoints #each point is a taylor
        push!(temp, sol.savedDers[k][i].coeffs[1])
    end
    display(plot!(p1,sol.savedTimes, temp,marker=(:circle),markersize=3,title="$(sol.algName)  _$(sol.absQ)_$(sol.simulStepCount)",label="x$k"#= ,xlims=(877.5,877.8),ylims=(-0.01,0.03) =#))
  end
    println("press enter to exit")
    readline() 
end =#
