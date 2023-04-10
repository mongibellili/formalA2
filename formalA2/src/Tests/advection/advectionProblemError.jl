using formalA2
#using OrdinaryDiffEq
#using ODEInterfaceDiffEq
using BSON
using StaticArrays

#using formalqssA
#using BenchmarkTools
#using Plots
#using OrdinaryDiffEq

#include("D://Advection.jl") 





#= function testN10()
    BSON.@load "solAdvection_N10_mliqss2e-6.bson" solmliqss2Interp
    BSON.@load "solVector_advection_N10_Rodas5Pe-9.bson" roads5VectorN10
    #err=getErrorByRodas(roads5VectorN10,solmliqss2Interp,1)
    #err=getAllErrorsByRodas(roads5VectorN10,solmliqss2Interp)
    err=getAverageErrorByRodas(roads5VectorN10,solmliqss2Interp)
    @show err

end =#

#testN10()

function testN100()
   # BSON.@load "solAdvection_N100_mliqss2e-6.bson" solmliqss2Interp
    BSON.@load "solAdvection_N100_NewLiqss-5.bson" solmliqss2Interp
   # BSON.@load "solAdvection_N100_mLiqssOlde-5.bson" solmliqss2Interp
  #  @show solmliqss2Interp.totalSteps #93498 #106703  #37807 ae flag index
   # BSON.@load "solVectAdvection_N100_Rodas5Pe-9.bson" roads5VectorN100
    BSON.@load "bson_base/solVectAdvection_N100_Feagin14e-12.bson" solFeagin14VectorN100
    #err=getErrorByRodas(roads5VectorN10,solmliqss2Interp,1)
    #err=getAllErrorsByRodas(roads5VectorN10,solmliqss2Interp)
    err=getAverageErrorByRodas(solFeagin14VectorN100,solmliqss2Interp)
    @show err  #err = ae0.0005370890318563965  # oldliqss2:err = 0.00048655207541598337  #  0.0005488222664032376 flag ae
    #for quan=1e-5  
        #err = 1.1492309021658802e-5  old: #1.0990442020733468e-5(elap) or err = 1.140279107860381e-5 AE###mliqss err = 1.8492865783211842e-5 

       #ujj=x- && removed abs(aii)>1e-6
       #AE: err = 1.5288030702973824e-5  :34608
       #AE flag: err = 1.0783420203099353e-5  :48766

     #time intgration=1.71s #2.1s
   #=   p1=getPlot(solmliqss2Interp,5)
     p1=plot!(p1,solFeagin14VectorN100[5],xlims=(0.0,3.5) ,ylims=(0.99,1.01)) =#
    # savefig(p1, "mliqss2AE_radau5_5")
   

     #= p2=getPlot(solmliqss2Interp,52)
  
     p2=plot!(p2,solFeagin14VectorN100[52],xlims=(5.0,10.0) ,ylims=(0.998,1.001))
     savefig(p2, "mliqssAE_feagin_52") =#
    # @show solmliqss2Interp.simulStepCount
end

testN100()