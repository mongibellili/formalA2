using formalA2
using XLSX
using Plots
#using BenchmarkTools

#include("C:/Users/belli/rel/formalA2/formalA2/test/typeH.jl")
include("./typeK.jl")
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
    p=getPlot!(solnmliqss2)
    p=plot!(x1)
    p=plot!(x2)
    savefig(p, "plot_$(prob.prname)_$(solType)_$(V)_$(ft).png")
   
    er1=getError(solnmliqss2Interp,1,x1)  
    er2=getError(solnmliqss2Interp,2,x2) 
    # timenmliqss=@belapsed QSS_Solve($odeprob,nmliqss2(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=100.0)
    resnmliqss= ("$(solnmliqss2.sysName)",(er1+er2)/2,solnmliqss2.totalSteps,solnmliqss2.simulStepCount,timenmliqss)
   # @show resnmliqss
end


function mainTest(solver::formalA2.QSSAlgorithm{solType, V})where {V,solType} 
    absTol=1e-7
    relTol=1e-4
    ft=2.0
    
   #=  funs=   [K_1,K_2,K_3,K_4,K_5,K_6,K_7,K_8,K_9,K_10,K_11,K_12,K_13,K_14,K_15,K_16,K_17,K_18,K_19,K_20,K_21,K_22,K_23,K_24,K_25,K_26,K_27,K_28,K_29,K_30,K_31,K_32,K_33,K_34,K_35,K_36,K_37,K_38,K_39,K_40,K_41,K_42,K_43,K_44,K_45,K_46,K_47,K_48,K_49,K_50,K_51,K_52,K_53,K_54,K_55,K_56,K_57,K_58,K_59,K_60,K_61,K_62,K_63,K_64,K_65,K_66,K_67,K_68,K_69,K_70,K_71,K_72,K_73,K_74,K_75,K_76,K_77,K_78,K_79,K_80,K_81,K_82,K_83,K_84,K_85,K_86,K_87,K_88,K_89,K_90,K_91,K_92,K_93,K_94,K_95,K_96,K_97,K_98,K_99,K_100,K_101,K_102,K_103,K_104,K_105,K_106,K_107,K_108,K_109,K_110,K_111,K_112,K_113,K_114,K_115,K_116,K_117,K_118,K_119,K_120,K_121,K_122,K_123,K_124,K_125,K_126,K_127,K_128,K_129,K_130,K_131,K_132,K_133,K_134,K_135,K_136,K_137,K_138,K_139,K_140,K_141,K_142,K_143,K_144,K_145,K_146,K_147,K_148,K_149,K_150,K_151,K_152,K_153,K_154,K_155,K_156,K_157,K_158,K_159,K_160,K_161,K_162,K_163,K_164,K_165,K_166,K_167,K_168,K_169,K_170,K_171,K_172,K_173,K_174,K_175,K_176,K_177,K_178,K_179,K_180,K_181,K_182,K_183,K_184,K_185,K_186,K_187,K_188,K_189,K_190,K_191,K_192,K_193,K_194,K_195,K_196,K_197,K_198,K_199,K_200,K_201,K_202,K_203,K_204,K_205,K_206,K_207,K_208,K_209,K_210,K_211,K_212,K_213,K_214,K_215,K_216,K_217,K_218,K_219,K_220,K_221,K_222,K_223,K_224,K_225,K_226,K_227,K_228,K_229,K_230,K_231,K_232,K_233,K_234,K_235,K_236,K_237,K_238,K_239,K_240,K_241,K_242,K_243,K_244,K_245,K_246,K_247,K_248,K_249,K_250,K_251,K_252,K_253,K_254,K_255,K_256,K_257,K_258,K_259,K_260,K_261,K_262,K_263,K_264,K_265,K_266,K_267,K_268,K_269,K_270,K_271,K_272,K_273,K_274,K_275,K_276,K_277,K_278,K_279,K_280,K_281,K_282,K_283,K_284,K_285,K_286,K_287,K_288,K_289,K_290,K_291,K_292,K_293,K_294,K_295,K_296,K_297,K_298,K_299,K_300,K_301,K_302,K_303,K_304,K_305,K_306,K_307,K_308,K_309,K_310,K_311,K_312,K_313,K_314,K_315,K_316,K_317,K_318,K_319,K_320,K_321,K_322,K_323,K_324,K_325,K_326,K_327,K_328,K_329,K_330,K_331,K_332,K_333,K_334,K_335,K_336,K_337,K_338,K_339,K_340,K_341,K_342,K_343,K_344,K_345,K_346,K_347,K_348,K_349,K_350,K_351,K_352,K_353,K_354,K_355,K_356,K_357,K_358,K_359,K_360,K_361,K_362,K_363,K_364,K_365,K_366,K_367,K_368,K_369,K_370,K_371,K_372,K_373,K_374,K_375,K_376,K_377,K_378,K_379,K_380,K_381,K_382,K_383,K_384,K_385,K_386,K_387,K_388,K_389,K_390,K_391,K_392,K_393,K_394,K_395,K_396,K_397,K_398,K_399,K_400,K_401,K_402,K_403,K_404,K_405,K_406,K_407,K_408,K_409,K_410,K_411,K_412,K_413,K_414,K_415,K_416,K_417,K_418,K_419,K_420,K_421,K_422,K_423,K_424,K_425,K_426,K_427,K_428,K_429,K_430,K_431,K_432,K_433,K_434,K_435,K_436,K_437,K_438,K_439,K_440,K_441,K_442,K_443,K_444,K_445,K_446,K_447,K_448,K_449,K_450,K_451,K_452,K_453,K_454,K_455,K_456,K_457,K_458,K_459,K_460,K_461,K_462,K_463,K_464,K_465,K_466,K_467,K_468,K_469,K_470,K_471,K_472,K_473,K_474,K_475,K_476,K_477,K_478,K_479,K_480,K_481,K_482,K_483,K_484,K_485,K_486,K_487,K_488,K_489,K_490,K_491,K_492,K_493,K_494,K_495,K_496,K_497,K_498,K_499,K_500,K_501,K_502,K_503,K_504,K_505,K_506,K_507,K_508,K_509,K_510,K_511,K_512,K_513,K_514,K_515,K_516,K_517,K_518,K_519,K_520,K_521,K_522,K_523,K_524,K_525,K_526,K_527,K_528,K_529,K_530,K_531,K_532,K_533,K_534,K_535,K_536,K_537,K_538,K_539,K_540,K_541,K_542,K_543,K_544,K_545,K_546,K_547,K_548,K_549,K_550,K_551,K_552,K_553,K_554,K_555,K_556,K_557,K_558,K_559,K_560,K_561,K_562,K_563,K_564,K_565,K_566,K_567,K_568,K_569,K_570,K_571,K_572,K_573,K_574,K_575,K_576,K_577,K_578,K_579,K_580,K_581,K_582,K_583,K_584,]
    =#
      funs=   [K_81]
   
    results=[]
    cntr=0
    for fun in funs
      cntr+=1
      if 100<cntr <102 @show fun end
      if 200<cntr <202 @show fun end
      if 300<cntr <302 @show fun end
      if 400<cntr <402 @show fun end
      if 500<cntr <502 @show fun end
      push!(results,solveProblem(fun,ft,solver,absTol,relTol))
    end
    @show results

   #=  XLSX.openxlsx("LTI_$(funs[1])_$(solType)_$V.xlsx", mode="w") do xf
      sheet = xf[1]
      sheet["A1"] = "LTI_H2"
      sheet["A2"] = "$(solType)_relTol=$relTol"
      sheet["A3"] = collect(("problem","error","totalSteps","simul_steps","time"))
      for i=4:length(funs)+3
        sheet["A$i"] = collect(results[i-3])
      end
    end =#
end


mainTest(nmliqss1())#normal analytic
#mainTest(nmliqss2())#normal analytic
#mainTest(mliqss1())#golden search analytic
#mainTest(nliqss1()) #iters
#mainTest(mliqssBounds1()) #analyticBounds