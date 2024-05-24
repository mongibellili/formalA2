using formalA2
using XLSX
#using BenchmarkTools

include("typeF.jl")
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

function mainTest1(solver1::formalA2.QSSAlgorithm{solType1, V},solver2::formalA2.QSSAlgorithm{solType2, V})where {V,solType1,solType2} 
    absTol=1e-6
    relTol=1e-3
    ft=100.0
    #funs=[F_1,F_2,F_3,F_4,F_5,F_6,F_7,F_8,F_9,F_10,F_11,F_12,F_13,F_14,F_15,F_16,F_17,F_18,F_19,F_20,F_21,F_22,F_23,F_24,F_25,F_26,F_27,F_28,F_29,F_30,F_31,F_32,F_33,F_34,F_35,F_36,F_37,F_38,F_39,F_40,F_41,F_42,F_43,F_44,F_45,F_46,F_47,F_48,F_49,F_50,F_51,F_52,F_53,F_54,F_55,F_56,F_57,F_58,F_59,F_60,]
    #funs=[F_1,F_2]
  
   funs=[F_776,F_777,F_778,F_779,F_780,F_781,F_782,F_783,F_784,F_785,F_786,F_787,F_788,F_789,F_790,F_791,F_792,F_793,F_794,F_795,F_796,F_797,F_798,F_799,F_800,F_801,F_802,F_803,F_804,F_805,F_806,F_807,F_808,F_809,F_810,F_811,F_812,F_813,F_814,F_815,F_816,F_817,F_818,F_819,F_820,F_821,F_822,F_823,F_824,F_825,F_826,F_827,F_828,F_829,F_830,F_831,F_832,F_833,F_834,F_835,F_836,F_837,F_838,F_839,F_840,F_841,F_842,F_843,F_844,F_845,F_846,F_847,F_848,F_849,F_850,F_851,F_852,F_853,F_854,F_855,F_856,F_857,F_858,F_859,F_860,F_861,F_862,F_863,F_864,F_865,F_866,F_867,F_868,F_869,F_870,F_871,F_872,F_873,F_874,F_875,F_876,F_877,F_878,F_879,F_880,F_881,F_882,F_883,F_884,F_885,F_886,F_887,F_888,F_889,F_890,F_891,F_892,F_893,F_894,F_895,F_896,F_897,F_898,F_899,F_900,F_901,F_902,F_903,F_904,F_905,F_906,F_907,F_908,F_909,F_910,F_911,F_912,F_913,F_914,F_915,F_916,F_917,F_918,F_919,F_920,F_921,F_922,F_923,F_924,F_925,F_926,F_927,F_928,F_929,F_930,F_931,F_932,F_933,F_934,F_935,F_936,F_937,F_938,F_939,F_940,F_941,F_942,F_943,F_944,F_945,F_946,F_947,F_948,F_949,F_950,F_951,F_952,F_953,F_954,F_955,F_956,F_957,F_958,F_959,F_960,F_961,F_962,F_963,F_964,F_965,F_966,F_967,F_968,F_969,F_970,F_971,F_972,F_973,F_974,F_975,F_976,F_977,F_978,F_979,F_980,F_981,F_982,F_983,F_984,F_985,F_986,F_987,F_988,F_989,F_990,F_991]
   
 
    results1=[];results2=[]
    for fun in funs
      sol1=solveProblem(fun,ft,solver1,absTol,relTol)
      sol2=solveProblem(fun,ft,solver2,absTol,relTol)
      if abs(sol1[3]-sol2[3])>3
        push!(results1,sol1)
        push!(results2,sol2)
      end
    end
    #@show results
    XLSX.openxlsx("LTI_F1_$(solType1)_$(solType2).xlsx", mode="w") do xf
      sheet = xf[1]
      sheet["A1"] = "LTI_F"
      sheet["A2"] = "$(solType1)_relTol=$relTol"
      sheet["A3"] = collect(("problem","error","totalSteps","simul_steps","time"))
      for i=4:length(results1)+3
        sheet["A$i"] = collect(results1[i-3])
      end

      #sheet["G1"] = "LTI_F"
      sheet["G2"] = "$(solType2)_relTol=$relTol"
      sheet["G3"] = collect(("problem","error","totalSteps","simul_steps","time"))
      for i=4:length(results2)+3
        sheet["G$i"] = collect(results2[i-3])
      end
    end
end

function mainTest2(solver1::formalA2.QSSAlgorithm{solType1, V},solver2::formalA2.QSSAlgorithm{solType2, V})where {V,solType1,solType2} 
  absTol=1e-6
  relTol=1e-3
  ft=100.0
  #funs=[F_1,F_2,F_3,F_4,F_5,F_6,F_7,F_8,F_9,F_10,F_11,F_12,F_13,F_14,F_15,F_16,F_17,F_18,F_19,F_20,F_21,F_22,F_23,F_24,F_25,F_26,F_27,F_28,F_29,F_30,F_31,F_32,F_33,F_34,F_35,F_36,F_37,F_38,F_39,F_40,F_41,F_42,F_43,F_44,F_45,F_46,F_47,F_48,F_49,F_50,F_51,F_52,F_53,F_54,F_55,F_56,F_57,F_58,F_59,F_60,]
  #funs=[F_1,F_2]

 
 funs=[F_992,F_993,F_994,F_995,F_996,F_997,F_998,F_999,F_1000,F_1001,F_1002,F_1003]
 #,F_1004,F_1005,F_1006,F_1007,F_1008,F_1009,F_1010,F_1011,F_1012,F_1013,F_1014,F_1015,F_1016,F_1017,F_1018,F_1019,F_1020,F_1021,F_1022,F_1023,F_1024,F_1025,F_1026,F_1027,F_1028,F_1029,F_1030,F_1031,F_1032,F_1033,F_1034,F_1035,F_1036,F_1037,F_1038,F_1039,F_1040,F_1041,F_1042,F_1043,F_1044,F_1045,F_1046,F_1047,F_1048,F_1049,F_1050,F_1051,F_1052,F_1053,F_1054,F_1055,F_1056,F_1057,F_1058,F_1059,F_1060,F_1061,F_1062,F_1063,F_1064,F_1065,F_1066,F_1067,F_1068,F_1069,F_1070,F_1071,F_1072,F_1073,F_1074,F_1075,F_1076,F_1077,F_1078,F_1079,F_1080,F_1081,F_1082,F_1083,F_1084,F_1085,F_1086,F_1087,F_1088,F_1089,F_1090,F_1091,F_1092,F_1093,F_1094,F_1095,F_1096]
 
  results1=[];results2=[]
  for fun in funs
    sol1=solveProblem(fun,ft,solver1,absTol,relTol)
    sol2=solveProblem(fun,ft,solver2,absTol,relTol)
    if abs(sol1[3]-sol2[3])>3
      push!(results1,sol1)
      push!(results2,sol2)
    end
  end
  #@show results
  XLSX.openxlsx("LTI_F2_$(solType1)_$(solType2).xlsx", mode="w") do xf
    sheet = xf[1]
    sheet["A1"] = "LTI_F"
    sheet["A2"] = "$(solType1)_relTol=$relTol"
    sheet["A3"] = collect(("problem","error","totalSteps","simul_steps","time"))
    for i=4:length(results1)+3
      sheet["A$i"] = collect(results1[i-3])
    end

    #sheet["G1"] = "LTI_F"
    sheet["G2"] = "$(solType2)_relTol=$relTol"
    sheet["G3"] = collect(("problem","error","totalSteps","simul_steps","time"))
    for i=4:length(results2)+3
      sheet["G$i"] = collect(results2[i-3])
    end
  end
end

function mainTest3(solver1::formalA2.QSSAlgorithm{solType1, V},solver2::formalA2.QSSAlgorithm{solType2, V})where {V,solType1,solType2} 
  absTol=1e-6
  relTol=1e-3
  ft=100.0
  #funs=[F_1,F_2,F_3,F_4,F_5,F_6,F_7,F_8,F_9,F_10,F_11,F_12,F_13,F_14,F_15,F_16,F_17,F_18,F_19,F_20,F_21,F_22,F_23,F_24,F_25,F_26,F_27,F_28,F_29,F_30,F_31,F_32,F_33,F_34,F_35,F_36,F_37,F_38,F_39,F_40,F_41,F_42,F_43,F_44,F_45,F_46,F_47,F_48,F_49,F_50,F_51,F_52,F_53,F_54,F_55,F_56,F_57,F_58,F_59,F_60,]
  #funs=[F_1,F_2]

 
 funs=[F_1097,F_1098,F_1099,F_1100,F_1101,F_1102,F_1103,F_1104,F_1105,F_1106,F_1107,F_1108,F_1109,F_1110,F_1111,F_1112,F_1113,F_1114,F_1115,F_1116,F_1117,F_1118,F_1119,F_1120,F_1121,F_1122,F_1123,F_1124,F_1125,F_1126,F_1127,F_1128,F_1129,F_1130,F_1131,F_1132,F_1133,F_1134,F_1135,F_1136,F_1137,F_1138,F_1139,F_1140,F_1141,F_1142,F_1143,F_1144,F_1145,F_1146,F_1147,F_1148,F_1149,F_1150,F_1151,F_1152,F_1153,F_1154,F_1155,F_1156,F_1157,F_1158,F_1159,F_1160,F_1161,F_1162,F_1163,F_1164,F_1165,F_1166,F_1167,F_1168,F_1169,F_1170,F_1171,F_1172,F_1173,F_1174,F_1175,F_1176,F_1177,F_1178,F_1179,F_1180,F_1181,F_1182,F_1183,F_1184,F_1185,F_1186,F_1187,F_1188,F_1189,F_1190,F_1191,F_1192,F_1193,F_1194,F_1195,F_1196,F_1197,F_1198,F_1199,F_1200,F_1201,F_1202,F_1203,F_1204,F_1205,F_1206,F_1207,F_1208,F_1209,F_1210,F_1211,F_1212,F_1213,F_1214,F_1215,F_1216,F_1217,F_1218,F_1219,F_1220,F_1221,F_1222,F_1223,F_1224,F_1225,F_1226,F_1227,F_1228,F_1229,F_1230,F_1231,F_1232,F_1233,F_1234,F_1235,F_1236,F_1237,F_1238,F_1239,F_1240,F_1241,F_1242,F_1243,F_1244,F_1245,F_1246,F_1247,F_1248,F_1249,F_1250,F_1251,F_1252,F_1253,F_1254,F_1255,F_1256,F_1257,F_1258,F_1259,F_1260,F_1261,F_1262,F_1263]
 
  results1=[];results2=[]
  for fun in funs
    sol1=solveProblem(fun,ft,solver1,absTol,relTol)
    sol2=solveProblem(fun,ft,solver2,absTol,relTol)
    if abs(sol1[3]-sol2[3])>3
      push!(results1,sol1)
      push!(results2,sol2)
    end
  end
  #@show results
  XLSX.openxlsx("LTI_F3_$(solType1)_$(solType2).xlsx", mode="w") do xf
    sheet = xf[1]
    sheet["A1"] = "LTI_F"
    sheet["A2"] = "$(solType1)_relTol=$relTol"
    sheet["A3"] = collect(("problem","error","totalSteps","simul_steps","time"))
    for i=4:length(results1)+3
      sheet["A$i"] = collect(results1[i-3])
    end

    #sheet["G1"] = "LTI_F"
    sheet["G2"] = "$(solType2)_relTol=$relTol"
    sheet["G3"] = collect(("problem","error","totalSteps","simul_steps","time"))
    for i=4:length(results2)+3
      sheet["G$i"] = collect(results2[i-3])
    end
  end
end

function mainTest4(solver1::formalA2.QSSAlgorithm{solType1, V},solver2::formalA2.QSSAlgorithm{solType2, V})where {V,solType1,solType2} 
  absTol=1e-6
  relTol=1e-3
  ft=100.0
  #funs=[F_1,F_2,F_3,F_4,F_5,F_6,F_7,F_8,F_9,F_10,F_11,F_12,F_13,F_14,F_15,F_16,F_17,F_18,F_19,F_20,F_21,F_22,F_23,F_24,F_25,F_26,F_27,F_28,F_29,F_30,F_31,F_32,F_33,F_34,F_35,F_36,F_37,F_38,F_39,F_40,F_41,F_42,F_43,F_44,F_45,F_46,F_47,F_48,F_49,F_50,F_51,F_52,F_53,F_54,F_55,F_56,F_57,F_58,F_59,F_60,]
  #funs=[F_1,F_2]

 
 funs=[F_1264,F_1265,F_1266,F_1267,F_1268,F_1269,F_1270,F_1271,F_1272,F_1273,F_1274,F_1275,F_1276,F_1277,F_1278,F_1279,F_1280,F_1281,F_1282,F_1283,F_1284,F_1285,F_1286,F_1287,F_1288,F_1289,F_1290,F_1291,F_1292,F_1293,F_1294,F_1295,F_1296,F_1297,F_1298,F_1299,F_1300,F_1301,F_1302,F_1303,F_1304,F_1305,F_1306,F_1307,F_1308,F_1309,F_1310,F_1311,F_1312,F_1313,F_1314,F_1315,F_1316,F_1317,F_1318,F_1319,F_1320,F_1321,F_1322,F_1323,F_1324,F_1325,F_1326,F_1327,F_1328,F_1329,F_1330,F_1331,F_1332,F_1333,F_1334,F_1335,F_1336,F_1337,F_1338,F_1339,F_1340,F_1341,F_1342,F_1343,F_1344,F_1345,F_1346,F_1347,F_1348,F_1349,F_1350,F_1351,F_1352,F_1353,F_1354,F_1355,F_1356,F_1357,F_1358,F_1359,F_1360,F_1361,F_1362,F_1363,F_1364,F_1365,F_1366,F_1367,F_1368,F_1369,F_1370,F_1371,F_1372,F_1373,F_1374,F_1375,F_1376,F_1377,F_1378,F_1379,F_1380,F_1381,F_1382,F_1383,F_1384,F_1385,F_1386,F_1387,F_1388,F_1389,F_1390,F_1391,F_1392,F_1393,F_1394,F_1395,F_1396,F_1397,F_1398,F_1399,F_1400,F_1401,F_1402,F_1403,F_1404,F_1405,F_1406,F_1407,F_1408,F_1409]
 
  results1=[];results2=[]
  for fun in funs
    sol1=solveProblem(fun,ft,solver1,absTol,relTol)
    sol2=solveProblem(fun,ft,solver2,absTol,relTol)
    if abs(sol1[3]-sol2[3])>3
      push!(results1,sol1)
      push!(results2,sol2)
    end
  end
  #@show results
  XLSX.openxlsx("LTI_F4_$(solType1)_$(solType2).xlsx", mode="w") do xf
    sheet = xf[1]
    sheet["A1"] = "LTI_F"
    sheet["A2"] = "$(solType1)_relTol=$relTol"
    sheet["A3"] = collect(("problem","error","totalSteps","simul_steps","time"))
    for i=4:length(results1)+3
      sheet["A$i"] = collect(results1[i-3])
    end

    #sheet["G1"] = "LTI_F"
    sheet["G2"] = "$(solType2)_relTol=$relTol"
    sheet["G3"] = collect(("problem","error","totalSteps","simul_steps","time"))
    for i=4:length(results2)+3
      sheet["G$i"] = collect(results2[i-3])
    end
  end
end

function mainTest5(solver1::formalA2.QSSAlgorithm{solType1, V},solver2::formalA2.QSSAlgorithm{solType2, V})where {V,solType1,solType2} 
  absTol=1e-6
  relTol=1e-3
  ft=100.0
  #funs=[F_1,F_2,F_3,F_4,F_5,F_6,F_7,F_8,F_9,F_10,F_11,F_12,F_13,F_14,F_15,F_16,F_17,F_18,F_19,F_20,F_21,F_22,F_23,F_24,F_25,F_26,F_27,F_28,F_29,F_30,F_31,F_32,F_33,F_34,F_35,F_36,F_37,F_38,F_39,F_40,F_41,F_42,F_43,F_44,F_45,F_46,F_47,F_48,F_49,F_50,F_51,F_52,F_53,F_54,F_55,F_56,F_57,F_58,F_59,F_60,]
  #funs=[F_1,F_2]

 
 funs=[F_1410,F_1411,F_1412,F_1413,F_1414,F_1415,F_1416,F_1417,F_1418,F_1419,F_1420,F_1421,F_1422,F_1423,F_1424,F_1425,F_1426,F_1427,F_1428,F_1429,F_1430,F_1431,F_1432,F_1433,F_1434,F_1435,F_1436,F_1437,F_1438,F_1439,F_1440,F_1441,F_1442,F_1443,F_1444,F_1445,F_1446,F_1447,F_1448,F_1449,F_1450,F_1451,F_1452,F_1453,F_1454,F_1455,F_1456,F_1457,F_1458,F_1459,F_1460,F_1461,F_1462,F_1463,F_1464,F_1465,F_1466,F_1467,F_1468,F_1469,F_1470,F_1471,F_1472,F_1473,F_1474,F_1475,F_1476,F_1477,F_1478,F_1479,F_1480,F_1481,F_1482,F_1483,F_1484,F_1485,F_1486,F_1487,F_1488,F_1489,F_1490,F_1491,F_1492,F_1493,F_1494,F_1495,F_1496,F_1497,F_1498,F_1499,F_1500,F_1501,F_1502,F_1503,F_1504,F_1505,F_1506,F_1507,F_1508,F_1509,F_1510,F_1511,F_1512,F_1513,F_1514,F_1515,F_1516,F_1517,F_1518,F_1519,F_1520,F_1521,F_1522,F_1523,F_1524,F_1525,F_1526,F_1527,F_1528,F_1529,F_1530,F_1531,F_1532,F_1533,F_1534,F_1535,F_1536,F_1537,F_1538,F_1539,F_1540,F_1541,F_1542,F_1543,F_1544,F_1545,F_1546,F_1547,F_1548,F_1549,F_1550,F_1551,F_1552,F_1553,F_1554,F_1555,F_1556,F_1557,F_1558,F_1559,F_1560,F_1561,F_1562,F_1563,F_1564,F_1565,F_1566,F_1567,F_1568,F_1569,F_1570,F_1571,F_1572,F_1573,F_1574,F_1575,F_1576]
 
  results1=[];results2=[]
  for fun in funs
    sol1=solveProblem(fun,ft,solver1,absTol,relTol)
    sol2=solveProblem(fun,ft,solver2,absTol,relTol)
    if abs(sol1[3]-sol2[3])>3
      push!(results1,sol1)
      push!(results2,sol2)
    end
  end
  #@show results
  XLSX.openxlsx("LTI_F5_$(solType1)_$(solType2).xlsx", mode="w") do xf
    sheet = xf[1]
    sheet["A1"] = "LTI_F"
    sheet["A2"] = "$(solType1)_relTol=$relTol"
    sheet["A3"] = collect(("problem","error","totalSteps","simul_steps","time"))
    for i=4:length(results1)+3
      sheet["A$i"] = collect(results1[i-3])
    end

    #sheet["G1"] = "LTI_F"
    sheet["G2"] = "$(solType2)_relTol=$relTol"
    sheet["G3"] = collect(("problem","error","totalSteps","simul_steps","time"))
    for i=4:length(results2)+3
      sheet["G$i"] = collect(results2[i-3])
    end
  end
end

function mainTest6(solver1::formalA2.QSSAlgorithm{solType1, V},solver2::formalA2.QSSAlgorithm{solType2, V})where {V,solType1,solType2} 
  absTol=1e-6
  relTol=1e-3
  ft=100.0
  #funs=[F_1,F_2,F_3,F_4,F_5,F_6,F_7,F_8,F_9,F_10,F_11,F_12,F_13,F_14,F_15,F_16,F_17,F_18,F_19,F_20,F_21,F_22,F_23,F_24,F_25,F_26,F_27,F_28,F_29,F_30,F_31,F_32,F_33,F_34,F_35,F_36,F_37,F_38,F_39,F_40,F_41,F_42,F_43,F_44,F_45,F_46,F_47,F_48,F_49,F_50,F_51,F_52,F_53,F_54,F_55,F_56,F_57,F_58,F_59,F_60,]
  #funs=[F_1,F_2]

 
 
 
 funs=[F_1723,F_1724,F_1725,F_1726,F_1727,F_1728,F_1729,F_1730,F_1731,F_1732,F_1733,F_1734,F_1735,F_1736,F_1737,F_1738,F_1739,F_1740,F_1741,F_1742,F_1743,F_1744,F_1745,F_1746,F_1747,F_1748,F_1749,F_1750,F_1751,F_1752,F_1753,F_1754,F_1755,F_1756,F_1757,F_1758,F_1759,F_1760,F_1761,F_1762,F_1763,F_1764,F_1765,F_1766,F_1767,F_1768,F_1769,F_1770,F_1771,F_1772,F_1773,F_1774,F_1775,F_1776,F_1777,F_1778,F_1779,F_1780,F_1781,F_1782,F_1783,F_1784,F_1785,F_1786,F_1787,F_1788,F_1789,F_1790,F_1791,F_1792,F_1793,F_1794,F_1795,F_1796,F_1797,F_1798,F_1799,F_1800,F_1801,F_1802,F_1803,F_1804,F_1805,F_1806,F_1807,F_1808,F_1809,F_1810,F_1811,F_1812,F_1813,F_1814,F_1815,F_1816,F_1817,F_1818,F_1819,F_1820,F_1821,F_1822,F_1823,F_1824,F_1825,F_1826,F_1827,F_1828,F_1829,F_1830,F_1831,F_1832,F_1833,F_1834,F_1835,F_1836,F_1837,F_1838,F_1839,F_1840,F_1841,F_1842,F_1843,F_1844,F_1845,F_1846,F_1847,F_1848,F_1849,F_1850,F_1851,F_1852,F_1853,F_1854,F_1855,F_1856,F_1857,F_1858,F_1859,F_1860,F_1861,F_1862,F_1863,F_1864,F_1865,F_1866,F_1867,F_1868,F_1869,F_1870,F_1871,F_1872,F_1873,F_1874,F_1875,F_1876,F_1877,F_1878,F_1879,F_1880,F_1881,F_1882,F_1883,F_1884,F_1885,F_1886,F_1887,F_1888,F_1889,F_1890,F_1891,F_1892,F_1893,F_1894,F_1895,F_1896,F_1897,F_1898,F_1899,F_1900,F_1901,F_1902,F_1903,F_1904,F_1905,F_1906,F_1907,F_1908,F_1909,F_1910,F_1911,F_1912,F_1913,F_1914,F_1915,F_1916,F_1917,F_1918,F_1919,F_1920,F_1921,F_1922,F_1923,F_1924,F_1925,F_1926,F_1927,F_1928,F_1929,F_1930,F_1931]
 

  results1=[];results2=[]
  for fun in funs
    sol1=solveProblem(fun,ft,solver1,absTol,relTol)
    sol2=solveProblem(fun,ft,solver2,absTol,relTol)
    if abs(sol1[3]-sol2[3])>3
      push!(results1,sol1)
      push!(results2,sol2)
    end
  end
  #@show results
  XLSX.openxlsx("LTI_F6_$(solType1)_$(solType2).xlsx", mode="w") do xf
    sheet = xf[1]
    sheet["A1"] = "LTI_F"
    sheet["A2"] = "$(solType1)_relTol=$relTol"
    sheet["A3"] = collect(("problem","error","totalSteps","simul_steps","time"))
    for i=4:length(results1)+3
      sheet["A$i"] = collect(results1[i-3])
    end

    #sheet["G1"] = "LTI_F"
    sheet["G2"] = "$(solType2)_relTol=$relTol"
    sheet["G3"] = collect(("problem","error","totalSteps","simul_steps","time"))
    for i=4:length(results2)+3
      sheet["G$i"] = collect(results2[i-3])
    end
  end
end

function mainTest7(solver1::formalA2.QSSAlgorithm{solType1, V},solver2::formalA2.QSSAlgorithm{solType2, V})where {V,solType1,solType2} 
  absTol=1e-6
  relTol=1e-3
  ft=100.0
  #funs=[F_1,F_2,F_3,F_4,F_5,F_6,F_7,F_8,F_9,F_10,F_11,F_12,F_13,F_14,F_15,F_16,F_17,F_18,F_19,F_20,F_21,F_22,F_23,F_24,F_25,F_26,F_27,F_28,F_29,F_30,F_31,F_32,F_33,F_34,F_35,F_36,F_37,F_38,F_39,F_40,F_41,F_42,F_43,F_44,F_45,F_46,F_47,F_48,F_49,F_50,F_51,F_52,F_53,F_54,F_55,F_56,F_57,F_58,F_59,F_60,]
  #funs=[F_1,F_2]

 
 
 
 funs=[F_1932,F_1933,F_1934,F_1935,F_1936,F_1937,F_1938,F_1939,F_1940,F_1941,F_1942,F_1943,F_1944,F_1945,F_1946,F_1947,F_1948,F_1949,F_1950,F_1951,F_1952,F_1953,F_1954,F_1955,F_1956,F_1957,F_1958,F_1959,F_1960,F_1961,F_1962,F_1963,F_1964,F_1965,F_1966,F_1967,F_1968,F_1969,F_1970,F_1971,F_1972,F_1973,F_1974,F_1975,F_1976,F_1977,F_1978,F_1979,F_1980,F_1981,F_1982,F_1983,F_1984,F_1985,F_1986,F_1987,F_1988,F_1989,F_1990,F_1991,F_1992,F_1993,F_1994,F_1995,F_1996,F_1997,F_1998,F_1999,F_2000,F_2001,F_2002,F_2003,F_2004,F_2005,F_2006,F_2007,F_2008,F_2009,F_2010,F_2011,F_2012,F_2013,F_2014,F_2015,F_2016,F_2017,F_2018,F_2019,F_2020,F_2021,F_2022,F_2023,F_2024,F_2025,F_2026,F_2027,F_2028,F_2029,F_2030,F_2031,F_2032,F_2033,F_2034,F_2035,F_2036,F_2037,F_2038,F_2039,F_2040,F_2041,F_2042,F_2043,F_2044,F_2045,F_2046,F_2047,F_2048,F_2049,F_2050,F_2051,F_2052,F_2053,F_2054,F_2055,F_2056,F_2057,F_2058,F_2059,F_2060,F_2061,F_2062,F_2063,F_2064,F_2065,F_2066,F_2067,F_2068,F_2069,F_2070,F_2071,F_2072,F_2073,F_2074,F_2075,F_2076,F_2077,F_2078,F_2079,F_2080,F_2081,F_2082,F_2083,F_2084,F_2085,F_2086,F_2087,F_2088,F_2089,F_2090,F_2091,F_2092,F_2093,F_2094,F_2095,F_2096,F_2097,F_2098,F_2099,F_2100,F_2101,F_2102,F_2103,F_2104,F_2105,F_2106,F_2107,F_2108,F_2109,F_2110,F_2111,F_2112,F_2113,F_2114,F_2115,F_2116,F_2117,F_2118,F_2119,F_2120,F_2121,F_2122,F_2123,F_2124,F_2125,F_2126,F_2127,F_2128,F_2129,F_2130,F_2131,F_2132,F_2133,F_2134,F_2135,F_2136,F_2137,F_2138,F_2139,F_2140,F_2141,F_2142,F_2143,F_2144,F_2145,F_2146,F_2147,F_2148,F_2149,F_2150,F_2151,F_2152,F_2153,F_2154,F_2155,F_2156,F_2157,F_2158,F_2159,F_2160,F_2161,F_2162,F_2163,F_2164,F_2165,F_2166,F_2167,F_2168,F_2169,F_2170,F_2171,F_2172,F_2173,F_2174,F_2175,F_2176,F_2177,F_2178,F_2179,F_2180,F_2181,F_2182]
 
 

  results1=[];results2=[]
  for fun in funs
    sol1=solveProblem(fun,ft,solver1,absTol,relTol)
    sol2=solveProblem(fun,ft,solver2,absTol,relTol)
    if abs(sol1[3]-sol2[3])>3
      push!(results1,sol1)
      push!(results2,sol2)
    end
  end
  #@show results
  XLSX.openxlsx("LTI_F7_$(solType1)_$(solType2).xlsx", mode="w") do xf
    sheet = xf[1]
    sheet["A1"] = "LTI_F"
    sheet["A2"] = "$(solType1)_relTol=$relTol"
    sheet["A3"] = collect(("problem","error","totalSteps","simul_steps","time"))
    for i=4:length(results1)+3
      sheet["A$i"] = collect(results1[i-3])
    end

    #sheet["G1"] = "LTI_F"
    sheet["G2"] = "$(solType2)_relTol=$relTol"
    sheet["G3"] = collect(("problem","error","totalSteps","simul_steps","time"))
    for i=4:length(results2)+3
      sheet["G$i"] = collect(results2[i-3])
    end
  end
end
function mainTest8(solver1::formalA2.QSSAlgorithm{solType1, V},solver2::formalA2.QSSAlgorithm{solType2, V})where {V,solType1,solType2} 
  absTol=1e-6
  relTol=1e-3
  ft=100.0
  #funs=[F_1,F_2,F_3,F_4,F_5,F_6,F_7,F_8,F_9,F_10,F_11,F_12,F_13,F_14,F_15,F_16,F_17,F_18,F_19,F_20,F_21,F_22,F_23,F_24,F_25,F_26,F_27,F_28,F_29,F_30,F_31,F_32,F_33,F_34,F_35,F_36,F_37,F_38,F_39,F_40,F_41,F_42,F_43,F_44,F_45,F_46,F_47,F_48,F_49,F_50,F_51,F_52,F_53,F_54,F_55,F_56,F_57,F_58,F_59,F_60,]
  #funs=[F_1,F_2]

 
 

 
 funs=[F_2183,F_2184,F_2185,F_2186,F_2187,F_2188,F_2189,F_2190,F_2191,F_2192,F_2193,F_2194,F_2195,F_2196,F_2197,F_2198,F_2199,F_2200,F_2201,F_2202,F_2203,F_2204,F_2205,F_2206,F_2207,F_2208,F_2209,F_2210,F_2211,F_2212,F_2213,F_2214,F_2215,F_2216,F_2217,F_2218,F_2219,F_2220,F_2221,F_2222,F_2223,F_2224,F_2225,F_2226,F_2227,F_2228,F_2229,F_2230,F_2231,F_2232,F_2233,F_2234,F_2235,F_2236,F_2237,F_2238,F_2239,F_2240,F_2241,F_2242,F_2243,F_2244,F_2245,F_2246,F_2247,F_2248,F_2249,F_2250,F_2251,F_2252,F_2253,F_2254,F_2255,F_2256,F_2257,F_2258,F_2259,F_2260,F_2261,F_2262,F_2263,F_2264,F_2265,F_2266,F_2267,F_2268,F_2269,F_2270,F_2271,F_2272,F_2273,F_2274,F_2275,F_2276,F_2277,F_2278,F_2279,F_2280,F_2281,F_2282,F_2283,F_2284,F_2285,F_2286,F_2287,F_2288,F_2289,F_2290,F_2291,F_2292,F_2293,F_2294,F_2295,F_2296,F_2297,F_2298,F_2299,F_2300,F_2301,F_2302,F_2303,F_2304,F_2305,F_2306,F_2307,F_2308,F_2309,F_2310,F_2311,F_2312,F_2313,F_2314,F_2315,F_2316,F_2317,F_2318,F_2319,F_2320,F_2321,F_2322,F_2323,F_2324,F_2325,F_2326,F_2327,F_2328,F_2329,F_2330,F_2331,F_2332,F_2333,F_2334,F_2335,F_2336,F_2337,F_2338,F_2339,F_2340,F_2341,F_2342,F_2343,F_2344,F_2345,F_2346,F_2347,F_2348,F_2349,F_2350,F_2351,F_2352,F_2353,F_2354,F_2355,F_2356,F_2357,F_2358,F_2359,F_2360,F_2361,F_2362,F_2363,F_2364,F_2365,F_2366,F_2367,F_2368,F_2369,F_2370,F_2371,F_2372,F_2373,F_2374,F_2375,F_2376,F_2377,F_2378,F_2379,F_2380,F_2381,F_2382,F_2383,F_2384,F_2385,F_2386,F_2387,F_2388,F_2389,F_2390,F_2391,F_2392,F_2393,F_2394,F_2395,F_2396,F_2397,F_2398,F_2399,F_2400,F_2401,F_2402,F_2403,F_2404,F_2405,F_2406,F_2407,F_2408,F_2409,F_2410,F_2411,F_2412,F_2413,F_2414,F_2415,F_2416,F_2417,F_2418,F_2419,F_2420,F_2421,F_2422,F_2423,F_2424,F_2425,F_2426,F_2427,F_2428,F_2429,F_2430,F_2431,F_2432,F_2433,F_2434,F_2435,F_2436,F_2437,F_2438,F_2439,F_2440,F_2441,F_2442,F_2443,F_2444,F_2445,F_2446,F_2447,F_2448,F_2449,F_2450,F_2451,F_2452,F_2453,F_2454]
 
 

  results1=[];results2=[]
  for fun in funs
    sol1=solveProblem(fun,ft,solver1,absTol,relTol)
    sol2=solveProblem(fun,ft,solver2,absTol,relTol)
    if abs(sol1[3]-sol2[3])>3
      push!(results1,sol1)
      push!(results2,sol2)
    end
  end
  #@show results
  XLSX.openxlsx("LTI_F8_$(solType1)_$(solType2).xlsx", mode="w") do xf
    sheet = xf[1]
    sheet["A1"] = "LTI_F"
    sheet["A2"] = "$(solType1)_relTol=$relTol"
    sheet["A3"] = collect(("problem","error","totalSteps","simul_steps","time"))
    for i=4:length(results1)+3
      sheet["A$i"] = collect(results1[i-3])
    end

    #sheet["G1"] = "LTI_F"
    sheet["G2"] = "$(solType2)_relTol=$relTol"
    sheet["G3"] = collect(("problem","error","totalSteps","simul_steps","time"))
    for i=4:length(results2)+3
      sheet["G$i"] = collect(results2[i-3])
    end
  end
end

function mainTest9(solver1::formalA2.QSSAlgorithm{solType1, V},solver2::formalA2.QSSAlgorithm{solType2, V})where {V,solType1,solType2} 
  absTol=1e-6
  relTol=1e-3
  ft=100.0
  #funs=[F_1,F_2,F_3,F_4,F_5,F_6,F_7,F_8,F_9,F_10,F_11,F_12,F_13,F_14,F_15,F_16,F_17,F_18,F_19,F_20,F_21,F_22,F_23,F_24,F_25,F_26,F_27,F_28,F_29,F_30,F_31,F_32,F_33,F_34,F_35,F_36,F_37,F_38,F_39,F_40,F_41,F_42,F_43,F_44,F_45,F_46,F_47,F_48,F_49,F_50,F_51,F_52,F_53,F_54,F_55,F_56,F_57,F_58,F_59,F_60,]
  #funs=[F_1,F_2]

 
 

 
 
 
 funs=[F_2455,F_2456,F_2457,F_2458,F_2459,F_2460,F_2461,F_2462,F_2463,F_2464,F_2465,F_2466,F_2467,F_2468,F_2469,F_2470,F_2471,F_2472,F_2473,F_2474,F_2475,F_2476,F_2477,F_2478,F_2479,F_2480,F_2481,F_2482,F_2483,F_2484,F_2485,F_2486,F_2487,F_2488,F_2489,F_2490,F_2491,F_2492,F_2493,F_2494,F_2495,F_2496,F_2497,F_2498,F_2499,F_2500,F_2501,F_2502,F_2503,F_2504,F_2505,F_2506,F_2507,F_2508,F_2509,F_2510,F_2511,F_2512,F_2513,F_2514,F_2515,F_2516,F_2517,F_2518,F_2519,F_2520,F_2521,F_2522,F_2523,F_2524,F_2525,F_2526,F_2527,F_2528,F_2529,F_2530,F_2531,F_2532,F_2533,F_2534,F_2535,F_2536,F_2537,F_2538,F_2539,F_2540,F_2541,F_2542,F_2543,F_2544,F_2545,F_2546,F_2547,F_2548,F_2549,F_2550,F_2551,F_2552,F_2553,F_2554,F_2555,F_2556,F_2557,F_2558,F_2559,F_2560,F_2561,F_2562,F_2563,F_2564,F_2565,F_2566,F_2567,F_2568,F_2569,F_2570,F_2571,F_2572,F_2573,F_2574,F_2575,F_2576,F_2577,F_2578,F_2579,F_2580,F_2581,F_2582,F_2583,F_2584,F_2585,F_2586,F_2587,F_2588,F_2589,F_2590,F_2591,F_2592,F_2593,F_2594,F_2595,F_2596,F_2597,F_2598,F_2599,F_2600,F_2601,F_2602,F_2603,F_2604,F_2605,F_2606,F_2607,F_2608,F_2609,F_2610,F_2611,F_2612,F_2613,F_2614,F_2615,F_2616,F_2617,F_2618,F_2619,F_2620,] 


  results1=[];results2=[]
  for fun in funs
    sol1=solveProblem(fun,ft,solver1,absTol,relTol)
    sol2=solveProblem(fun,ft,solver2,absTol,relTol)
    if abs(sol1[3]-sol2[3])>3
      push!(results1,sol1)
      push!(results2,sol2)
    end
  end
  #@show results
  XLSX.openxlsx("LTI_F9_$(solType1)_$(solType2).xlsx", mode="w") do xf
    sheet = xf[1]
    sheet["A1"] = "LTI_F"
    sheet["A2"] = "$(solType1)_relTol=$relTol"
    sheet["A3"] = collect(("problem","error","totalSteps","simul_steps","time"))
    for i=4:length(results1)+3
      sheet["A$i"] = collect(results1[i-3])
    end

    #sheet["G1"] = "LTI_F"
    sheet["G2"] = "$(solType2)_relTol=$relTol"
    sheet["G3"] = collect(("problem","error","totalSteps","simul_steps","time"))
    for i=4:length(results2)+3
      sheet["G$i"] = collect(results2[i-3])
    end
  end
end
#mainTest(nmliqss1())#normal analytic
#mainTest(mliqss1())#golden search analytic
#mainTest(nliqss1()) #iters
#= mainTest1(mliqssBounds1(),nmliqss1()) #analyticBounds
mainTest2(mliqssBounds1(),nmliqss1()) #analyticBounds
mainTest3(mliqssBounds1(),nmliqss1()) #analyticBounds
mainTest4(mliqssBounds1(),nmliqss1()) #analyticBounds
mainTest5(mliqssBounds1(),nmliqss1()) #analyticBounds
mainTest6(mliqssBounds1(),nmliqss1()) #analyticBounds
mainTest7(mliqssBounds1(),nmliqss1()) #analyticBounds
mainTest8(mliqssBounds1(),nmliqss1()) #analyticBounds
mainTest9(mliqssBounds1(),nmliqss1()) #analyticBounds =#

mainTest2(mliqss1(),nmliqss1()) #analytic_Optim