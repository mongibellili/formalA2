using formalA2
using XLSX
using Plots
#using BenchmarkTools

#include("C:/Users/belli/rel/formalA2/formalA2/test/typeH.jl")
include("./typeJ.jl")
function solveProblem(prblem::Function,ft::Float64,solver::formalA2.QSSAlgorithm{solType, V},absTol,relTol)where {V,solType} 
    pr=prblem()
    prob=pr[1]
    x1=pr[2]
    x2=pr[3]
    timenmliqss=0.0
    #= absTol=1e-5
    relTol=1e-3 =#
    tspan=(0.0,ft)
    solnmliqss2=solve(prob,solver,abstol=absTol,saveat=0.01,reltol=relTol,tspan)
    solnmliqss2Interp=solInterpolated(solnmliqss2,0.01)
   #=  p=getPlot!(solnmliqss2)
    p=plot!(x1)
    p=plot!(x2)
    savefig(p, "plotInterpl_$(solType)_$V.png") =#
   
    er1=getError(solnmliqss2Interp,1,x1)  
    er2=getError(solnmliqss2Interp,2,x2) 
    # timenmliqss=@belapsed QSS_Solve($odeprob,nmliqss2(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=100.0)
    resnmliqss= ("$(solnmliqss2.sysName)",(er1+er2)/2,solnmliqss2.totalSteps,solnmliqss2.simulStepCount,timenmliqss)
   # @show resnmliqss
end


function mainTest(solver::formalA2.QSSAlgorithm{solType, V})where {V,solType} 
    absTol=1e-5
    relTol=1e-3
    ft=3.0
    
  funs=[J_1,J_2,J_3,J_4,J_5,J_6,J_7,J_8,J_9,J_10,J_11,J_12,J_13,J_14,J_15,J_16,J_17,J_18,J_19,J_20,J_21,J_22,J_23,J_24,J_25,J_26,J_27,J_28,J_29,J_30,J_31,J_32,J_33,J_34,J_35,J_36,J_37,J_38,J_39,J_40,J_41,J_42,J_43,J_44,J_45,J_46,J_47,J_48,J_49,J_50,J_51,J_52,J_53,J_54,J_55,J_56,J_57,J_58,J_59,J_60,J_61,J_62,J_63,J_64,J_65,J_66,J_67,J_68,J_69,J_70,J_71,J_72,J_73,J_74,J_75,J_76,J_77,J_78,J_79,J_80,J_81,J_82,J_83,J_84,J_85,J_86,J_87,J_88,J_89,J_90,J_91,J_92,J_93,J_94,J_95,J_96,J_97,J_98,J_99,J_100,J_101,J_102,J_103,J_104,J_105,J_106,J_107,J_108,J_109,J_110,J_111,J_112,J_113,J_114,J_115,J_116,J_117,J_118,J_119,J_120,J_121,J_122,J_123,J_124,J_125,J_126,J_127,J_128,J_129,J_130,J_131,J_132,J_133,J_134,J_135,J_136,J_137,J_138,J_139,J_140,J_141,J_142,J_143,J_144,J_145,J_146,J_147,J_148,J_149,J_150,J_151,J_152,J_153,J_154,J_155,J_156,J_157,J_158,J_159,J_160,J_161,J_162,J_163,J_164,J_165,J_166,J_167,J_168,J_169,J_170,J_171,J_172,J_173,J_174,J_175,J_176,J_177,J_178,J_179,J_180,J_181,J_182,J_183,J_184,J_185,J_186,J_187,J_188,J_189,J_190,J_191,J_192,J_193,J_194,J_195,J_196,J_197,J_198,J_199,J_200,J_201,J_202,J_203,J_204,J_205,J_206,J_207,J_208,J_209,J_210,J_211,J_212,J_213,J_214,J_215,J_216,J_217,J_218,J_219,J_220,J_221,J_222,J_223,J_224,J_225,J_226,J_227,J_228,J_229,J_230,J_231,J_232,J_233,J_234,J_235,J_236,J_237,J_238,J_239,J_240,J_241,J_242,J_243,J_244,J_245,J_246,J_247,J_248,J_249,J_250,J_251,J_252,J_253,J_254,J_255,J_256,J_257,J_258,J_259,J_260,J_261,J_262,J_263,J_264,J_265,J_266,J_267,J_268,J_269,J_270,J_271,J_272,J_273,J_274,J_275,J_276,J_277,J_278,J_279,J_280,J_281,J_282,J_283,J_284,J_285,J_286,J_287,J_288,J_289,J_290,J_291,J_292,J_293,J_294,J_295,J_296,J_297,J_298,J_299,J_300,J_301,J_302,J_303,J_304,J_305,J_306,J_307,J_308,J_309,J_310,J_311,J_312,J_313,J_314,J_315,J_316,J_317,J_318,J_319,J_320,J_321,J_322,J_323,J_324,J_325,J_326,J_327,J_328,J_329,J_330,J_331,J_332,J_333,J_334,J_335,J_336,J_337,J_338,J_339,J_340,J_341,J_342,J_343,J_344,J_345,J_346,J_347,J_348,J_349,J_350,J_351,J_352,J_353,J_354,J_355,J_356,J_357,J_358,J_359,J_360,J_361,J_362,J_363,J_364,J_365,J_366,J_367,J_368,J_369,J_370,J_371,J_372,J_373,J_374,J_375,J_376,J_377,J_378,J_379,J_380,J_381,J_382,J_383,J_384,J_385,J_386,J_387,J_388,J_389,J_390,J_391,J_392,J_393,J_394,J_395,J_396,J_397,J_398,J_399,J_400,J_401,J_402,J_403,J_404,J_405,J_406,J_407,J_408,J_409,J_410,J_411,J_412,J_413,J_414,J_415,J_416,J_417,J_418,J_419,J_420,J_421,J_422,J_423,J_424,J_425,J_426,J_427,J_428,J_429,J_430,J_431,J_432,J_433,J_434,J_435,J_436,J_437,J_438,J_439,J_440,J_441,J_442,J_443,J_444,J_445,J_446,J_447,J_448,J_449,J_450,J_451,J_452,J_453,J_454,J_455,J_456,J_457,J_458,J_459,J_460,J_461,J_462,J_463,J_464,J_465,J_466,J_467,J_468,J_469,J_470,J_471,J_472,J_473,J_474,J_475,J_476,J_477,J_478,J_479,J_480,J_481,J_482,J_483,J_484,J_485,J_486,J_487,J_488,J_489,J_490,J_491,J_492,J_493,J_494,J_495,J_496,J_497,J_498,J_499,J_500,J_501,J_502,J_503,J_504,J_505,J_506,J_507,J_508,J_509,J_510,]

    results=[]
    cntr=0
    for fun in funs
      cntr+=1
      if 100<cntr <102 @show fun end
      if 200<cntr <202 @show fun end
      if 300<cntr <302 @show fun end
      if 400<cntr <402 @show fun end
      push!(results,solveProblem(fun,ft,solver,absTol,relTol))
    end
    #@show results
   # @show results
    XLSX.openxlsx("LTI_$(funs[1])_$(solType)_$V.xlsx", mode="w") do xf
      sheet = xf[1]
      sheet["A1"] = "LTI_H2"
      sheet["A2"] = "$(solType)_relTol=$relTol"
      sheet["A3"] = collect(("problem","error","totalSteps","simul_steps","time"))
      for i=4:length(funs)+3
        sheet["A$i"] = collect(results[i-3])
      end
    end
end


mainTest(nmliqss1())#normal analytic
#mainTest(nmliqss2())#normal analytic
#mainTest(mliqss1())#golden search analytic
#mainTest(nliqss1()) #iters
#mainTest(mliqssBounds1()) #analyticBounds