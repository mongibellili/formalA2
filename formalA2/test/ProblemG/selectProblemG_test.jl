using formalA2
using XLSX
#using BenchmarkTools
include("typeG.jl")
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
   # save_Sol(solnmliqss2)
    solnmliqss2Interp=solInterpolated(solnmliqss2,0.01)
    er1=getError(solnmliqss2Interp,1,x1)  
    er2=getError(solnmliqss2Interp,2,x2) 
    # timenmliqss=@belapsed QSS_Solve($odeprob,nmliqss2(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=100.0)
    resnmliqss= ("$(solnmliqss2.sysName)",(er1+er2)/2,solnmliqss2.totalSteps,solnmliqss2.simulStepCount,timenmliqss)
   # @show resnmliqss
end


function mainTest(solver::formalA2.QSSAlgorithm{solType, V})where {V,solType} 
    absTol=1e-6
    relTol=1e-3
    ft=10.0
    
   #funs=[G_1,G_2,G_3,G_4,G_5,G_6,G_7,G_8,G_9,G_10,G_11,G_12,G_13,G_14,G_15,G_16,G_17,G_18,G_19,G_20,G_21,G_22,G_23,G_24,G_25,G_26,G_27,G_28,G_29,G_30,G_31,G_32,G_33,G_34,G_35,G_36,G_37,G_38,G_39,G_40,G_41,G_42,G_43,G_44,G_45,G_46,G_47,G_48,G_49,G_50,G_51,G_52,G_53,G_54,G_55,G_56,G_57,G_58,G_59,G_60,G_61,G_62,G_63,G_64,G_65,G_66,G_67,G_68,G_69,G_70,G_71,G_72,G_73,G_74,G_75,G_76,G_77,G_78,G_79,G_80,G_81,G_82,G_83,G_84,G_85,G_86,G_87,G_88,G_89,G_90,G_91,G_92,G_93,G_94,G_95,G_96,G_97,G_98,G_99,G_100,G_101,G_102,G_103,G_104,G_105,G_106,G_107,G_108,G_109,G_110,G_111,G_112,G_113,G_114,G_115,G_116,G_117,G_118,G_119,G_120,G_121,G_122,G_123,G_124,G_125,G_126,G_127,G_128,G_129,G_130,G_131,G_132,G_133,G_134,G_135,G_136,G_137,G_138,G_139,G_140,G_141,G_142,G_143,G_144,G_145,G_146,G_147,G_148,G_149,G_150,G_151,G_152,G_153,G_154,G_155,G_156,G_157,G_158,G_159,G_160,G_161,G_162,G_163,G_164,G_165,G_166,G_167,G_168,G_169,G_170,G_171,G_172,G_173,G_174,G_175,G_176,G_177,G_178,G_179,G_180,G_181,G_182,G_183,G_184,G_185,G_186,G_187,G_188,G_189,G_190,G_191,G_192,G_193,G_194,G_195,G_196,G_197,G_198,G_199,G_200,G_201,G_202,G_203,G_204,G_205,G_206,G_207,G_208,G_209,G_210,G_211,G_212,G_213,G_214,G_215,G_216,G_217,G_218,G_219,G_220,G_221,G_222,G_223,G_224,G_225,G_226,G_227,G_228,G_229,G_230,G_231,G_232,G_233,G_234,G_235,G_236,G_237,G_238,G_239,G_240,G_241,G_242,G_243,G_244,G_245,G_246,G_247,G_248,G_249,G_250,G_251,G_252,G_253,G_254,G_255,G_256,G_257,G_258,G_259,G_260,G_261,G_262,G_263,G_264,G_265,G_266,G_267,G_268,G_269,G_270,G_271,G_272,G_273,G_274,G_275,G_276,G_277,G_278,G_279,G_280,G_281,G_282,G_283,G_284,G_285,G_286,G_287,G_288,G_289,G_290,G_291,G_292,G_293,G_294,G_295,G_296,G_297,G_298,G_299,G_300,G_301,G_302,G_303,G_304,G_305,G_306,G_307,G_308,G_309,G_310,G_311,G_312,G_313,G_314,G_315,G_316,G_317,G_318,G_319,G_320,G_321,G_322,G_323,G_324,G_325,G_326,G_327,G_328,G_329,G_330,G_331,G_332,G_333,G_334,G_335,G_336,G_337,G_338,G_339,G_340,G_341,G_342,G_343,G_344,G_345,G_346,G_347,G_348,G_349,G_350,G_351,G_352,G_353,G_354,G_355,G_356,G_357,G_358,G_359,G_360,G_361,G_362,G_363,G_364,G_365,G_366,G_367,G_368,G_369,G_370,G_371,G_372,G_373,G_374,G_375,G_376,G_377,G_378,G_379,G_380,G_381,G_382,G_383,G_384,G_385,G_386,G_387,G_388,G_389,G_390,G_391,G_392,G_393,G_394,G_395,G_396,G_397,G_398,G_399,G_400,G_401,G_402,G_403,G_404,G_405,G_406,G_407,G_408,G_409,G_410,G_411,G_412,G_413,G_414,G_415,G_416,G_417,G_418,G_419,G_420,G_421,G_422,G_423,G_424,G_425,G_426,G_427,G_428,G_429,G_430,G_431,G_432,G_433,G_434,G_435,G_436,G_437,G_438,G_439,G_440,G_441,G_442,G_443,G_444,G_445,G_446,G_447,G_448,G_449,G_450,G_451,G_452,G_453,G_454,G_455,G_456,G_457,G_458,G_459,G_460,G_461,G_462,G_463,G_464,G_465,G_466,G_467,G_468,G_469,G_470,G_471,G_472,G_473,G_474,G_475,G_476,G_477,G_478,G_479,G_480,G_481,G_482,G_483,G_484,G_485,G_486,G_487,G_488,G_489,G_490,G_491,G_492,G_493,G_494,G_495,G_496,G_497,G_498,G_499,G_500,G_501,G_502,G_503,G_504,G_505,G_506,G_507,G_508,G_509,G_510,G_511,G_512,G_513,G_514,G_515,G_516,G_517,G_518,G_519,G_520,G_521,G_522,G_523,G_524,G_525,G_526,G_527,G_528,G_529,G_530,G_531,G_532,G_533,G_534,G_535,G_536,G_537,G_538,G_539,G_540,G_541,G_542,G_543,G_544,G_545,G_546,G_547,G_548,G_549,G_550,G_551,G_552,G_553,G_554,G_555,G_556,G_557,G_558,G_559,G_560,G_561,G_562,G_563,G_564,G_565,G_566,G_567,G_568,G_569,G_570,G_571,G_572,G_573,G_574,G_575,G_576,G_577,G_578,G_579,G_580,G_581,G_582,G_583,G_584,G_585,G_586,G_587,G_588,G_589,G_590,G_591,G_592,G_593,G_594,G_595,G_596,G_597,G_598,G_599,G_600,G_601,G_602,G_603,G_604,G_605,G_606,G_607,G_608,G_609,G_610,G_611,G_612,G_613,G_614,G_615,G_616,G_617,G_618,G_619,G_620,G_621,G_622,G_623,G_624,G_625,G_626,G_627,G_628,G_629,G_630,G_631,G_632,G_633,G_634,G_635,G_636,G_637,G_638,G_639,G_640,G_641,G_642,G_643,G_644,G_645,G_646,G_647,G_648,G_649,G_650,G_651,G_652,G_653,G_654,G_655,G_656,G_657,G_658,G_659,G_660,G_661,G_662,G_663,G_664,G_665,G_666,G_667,G_668,G_669,G_670,G_671,G_672,G_673,G_674,G_675,G_676,G_677,G_678,G_679,G_680,G_681,G_682,G_683,G_684,G_685,G_686,G_687,G_688,G_689,G_690,G_691,G_692,G_693,G_694,G_695,G_696,G_697,G_698,G_699,G_700,G_701,G_702,G_703,G_704,G_705,G_706,G_707,G_708,G_709,G_710,]

   funs=[G_13]
    results=[]
    for fun in funs
      @show fun
      push!(results,solveProblem(fun,ft,solver,absTol,relTol))
    end
    @show results
  #=   @show results
    XLSX.openxlsx("LTI_GG_$(solType).xlsx", mode="w") do xf
      sheet = xf[1]
      sheet["A1"] = "LTI_GG"
      sheet["A2"] = "$(solType)_relTol=$relTol"
      sheet["A3"] = collect(("problem","error","totalSteps","simul_steps","time"))
      for i=4:length(funs)+3
        sheet["A$i"] = collect(results[i-3])
      end
    end =#
end


mainTest(nmliqss1())#normal analytic
#mainTest(mliqss1())#golden search analytic
#mainTest(nliqss1()) #iters
#mainTest(mliqssBounds1()) #analyticBounds