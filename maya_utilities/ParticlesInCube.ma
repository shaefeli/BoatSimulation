//Maya ASCII 2017 scene
//Name: ParticlesInCube.ma
//Last modified: Mon, Nov 27, 2017 10:22:21 AM
//Codeset: 1252
requires maya "2017";
currentUnit -l centimeter -a degree -t film;
fileInfo "application" "maya";
fileInfo "product" "Maya 2017";
fileInfo "version" "2017";
fileInfo "cutIdentifier" "201606150345-997974";
fileInfo "osv" "Microsoft Windows 7 Business Edition, 64-bit Windows 7 Service Pack 1 (Build 7601)\n";
fileInfo "license" "student";
createNode transform -s -n "persp";
	rename -uid "7185C08C-45B0-9914-5E7F-F58B088B941F";
	setAttr ".v" no;
	setAttr ".t" -type "double3" 28 21 28 ;
	setAttr ".r" -type "double3" -27.938352729602379 44.999999999999972 -5.172681101354183e-014 ;
createNode camera -s -n "perspShape" -p "persp";
	rename -uid "B7B76270-45C3-73B5-C146-C9945484BF41";
	setAttr -k off ".v" no;
	setAttr ".fl" 34.999999999999993;
	setAttr ".coi" 44.82186966202994;
	setAttr ".imn" -type "string" "persp";
	setAttr ".den" -type "string" "persp_depth";
	setAttr ".man" -type "string" "persp_mask";
	setAttr ".hc" -type "string" "viewSet -p %camera";
	setAttr ".ai_translator" -type "string" "perspective";
createNode transform -s -n "top";
	rename -uid "AD93D0F2-479F-6FCC-C4A5-0582C3675BBD";
	setAttr ".v" no;
	setAttr ".t" -type "double3" 0 1000.1 0 ;
	setAttr ".r" -type "double3" -89.999999999999986 0 0 ;
createNode camera -s -n "topShape" -p "top";
	rename -uid "AAFABFA7-4695-1258-6C75-8F93F3B1289B";
	setAttr -k off ".v" no;
	setAttr ".rnd" no;
	setAttr ".coi" 1000.1;
	setAttr ".ow" 30;
	setAttr ".imn" -type "string" "top";
	setAttr ".den" -type "string" "top_depth";
	setAttr ".man" -type "string" "top_mask";
	setAttr ".hc" -type "string" "viewSet -t %camera";
	setAttr ".o" yes;
	setAttr ".ai_translator" -type "string" "orthographic";
createNode transform -s -n "front";
	rename -uid "07423675-4AD0-D930-EFB8-D08CCD033B04";
	setAttr ".v" no;
	setAttr ".t" -type "double3" 0 0 1000.1 ;
createNode camera -s -n "frontShape" -p "front";
	rename -uid "9A7FC0F4-4EAD-EFA5-AD34-52A14EB2E365";
	setAttr -k off ".v" no;
	setAttr ".rnd" no;
	setAttr ".coi" 1000.1;
	setAttr ".ow" 30;
	setAttr ".imn" -type "string" "front";
	setAttr ".den" -type "string" "front_depth";
	setAttr ".man" -type "string" "front_mask";
	setAttr ".hc" -type "string" "viewSet -f %camera";
	setAttr ".o" yes;
	setAttr ".ai_translator" -type "string" "orthographic";
createNode transform -s -n "side";
	rename -uid "505A48E9-4E7B-864B-F0E8-BDA4AB8E4EE8";
	setAttr ".v" no;
	setAttr ".t" -type "double3" 1000.1 0 0 ;
	setAttr ".r" -type "double3" 0 89.999999999999986 0 ;
createNode camera -s -n "sideShape" -p "side";
	rename -uid "D11E21BD-4D87-E472-ECBC-819FC310152D";
	setAttr -k off ".v" no;
	setAttr ".rnd" no;
	setAttr ".coi" 1000.1;
	setAttr ".ow" 30;
	setAttr ".imn" -type "string" "side";
	setAttr ".den" -type "string" "side_depth";
	setAttr ".man" -type "string" "side_mask";
	setAttr ".hc" -type "string" "viewSet -s %camera";
	setAttr ".o" yes;
	setAttr ".ai_translator" -type "string" "orthographic";
createNode transform -n "pCube1";
	rename -uid "649B1D52-4142-B461-C245-F184D43846BB";
	setAttr ".s" -type "double3" 10 2 10 ;
createNode mesh -n "pCubeShape1" -p "pCube1";
	rename -uid "60008B7B-4CAB-75B2-C86E-43A99298F542";
	setAttr -k off ".v";
	setAttr ".vir" yes;
	setAttr ".vif" yes;
	setAttr ".pv" -type "double2" 0.5 0.375 ;
	setAttr ".uvst[0].uvsn" -type "string" "map1";
	setAttr ".cuvs" -type "string" "map1";
	setAttr ".dcc" -type "string" "Ambient+Diffuse";
	setAttr ".covm[0]"  0 1 1;
	setAttr ".cdvm[0]"  0 1 1;
	setAttr ".ai_translator" -type "string" "polymesh";
createNode nucleus -n "nucleus1";
	rename -uid "C11D1C29-44B6-A668-506C-24B574D4BE7C";
createNode transform -n "nParticle1";
	rename -uid "287A3003-4FDB-C9D2-628D-23864FCF4CEC";
createNode nParticle -n "nParticleShape1" -p "nParticle1";
	rename -uid "B44B925C-457C-8339-4764-FBA63FC2DCFD";
	addAttr -s false -ci true -sn "lifespanPP" -ln "lifespanPP" -dt "doubleArray";
	addAttr -ci true -h true -sn "lifespanPP0" -ln "lifespanPP0" -dt "doubleArray";
	addAttr -ci true -sn "lifespan" -ln "lifespan" -at "double";
	addAttr -s false -ci true -sn "rgbPP" -ln "rgbPP" -dt "vectorArray";
	addAttr -ci true -h true -sn "rgbPP0" -ln "rgbPP0" -dt "vectorArray";
	addAttr -s false -ci true -sn "opacityPP" -ln "opacityPP" -dt "doubleArray";
	addAttr -ci true -h true -sn "opacityPP0" -ln "opacityPP0" -dt "doubleArray";
	addAttr -s false -ci true -sn "radiusPP" -ln "radiusPP" -dt "doubleArray";
	addAttr -ci true -h true -sn "radiusPP0" -ln "radiusPP0" -dt "doubleArray";
	setAttr -k off ".v";
	setAttr ".gf" -type "Int32Array" 0 ;
	setAttr ".pos0" -type "vectorArray" 100 -0.66679660887238512 0 -0.71474334021684527 -1.2360656420273748
		 -1.1265540777550371 -0.016937259571844887 -0.99632859203697921 0.27370098797482573
		 -0.88794561732906696 0.53841723174491651 -0.2038082140009351 -1.4436959113275272 -0.76658828733700302
		 -0.20149885341577375 0.34965870688435174 -0.70095093607547976 -0.37581172530022972
		 -0.3874539296951276 -0.020862331261456979 -1.0563568635398497 -0.4939537798536186 -1.2680292229869103
		 1.4644903664312778 -0.85153612408773616 -0.0017703270955613526 0.48237712984882819
		 -0.61830779416651249 -1.3390400333773749 0.92835760753604291 -1.7520603255932521 -0.6684086658466406
		 0.91458061267919677 -1.7716812654654157 -1.4184234818881685 -0.52813601832464208
		 -0.3519437636236013 0.068315042278312199 -0.83205486576202936 -1.7766449822704977 -1.2380388539258824
		 -0.75530293228794676 -0.98367787008476526 0.057457458560629338 0.21989617779170773
		 -1.8482025950701464 -0.49306116681440737 0.069428625555718784 0.16207646442674895 -1.808226492874865
		 0.31498522269708362 0.067067229809000484 -0.17579139818696415 0.57923056777398763
		 -0.50429046501603025 -0.29368100198296854 0.049781302072580053 -2.1380021862317826 -0.68412025487131634
		 0.55517092975912308 -0.52689954938124761 -1.2541353580853687 0.13634228882056618
		 -1.1491828738781158 -1.9107260122711232 -0.63031270794121497 -0.31292506502572121 -1.3412564060910477
		 0.34787917164119908 -1.9731282431705985 0.51093631915196691 0.82111425888116973 -0.88210098481112387 -0.60528748626016338
		 -0.68614334831159463 0.25077707283018857 -0.37559023034335653 -0.22698788806660561
		 -1.4369542395907935 0.17562397325214152 0.77314639990084977 -0.64035979715701852 -0.41989169822119193
		 0.30530585487688316 -1.532897365682468 -0.65079106256056252 -0.39170126356229962
		 -0.62386646104324439 -0.73778021919335335 -0.3249822483705373 -1.7084702486676782 -1.2584171656659286
		 0.47124167277000312 -2.0967268991139951 -1.8705029381669016 -0.64599111222524364
		 -1.110223991968295 0.51625801278447692 0.84056474364042522 -1.2258002781417852 -0.72118791300002205
		 -0.29868958030343151 -1.8748246262446437 -0.4110363073098135 0.88318817315952181
		 -0.2666443813637876 0.74044723303497983 0.68414314712191526 -0.86783177389180322 -0.60017204746458153
		 0.81003502822377416 0.55210462240728697 0.16276905390473717 1.3675614288406166 -0.70336307173853707 -1.6593507784041599
		 -1.1614287611081238 -0.5330550658213532 -0.25731017226095715 0.87186712364927532
		 -1.5223784266345923 -0.7336457918427236 0.72703695554867531 -0.046697567402410489 -1.1600193343506184
		 -0.61452090730472841 -1.3911315067547811 -0.84675698441293934 0.52231099033010087
		 0.051864991010059924 -0.29129440325463563 -0.073131171764651989 -0.13142030252899206 -0.49550588729509781
		 -0.044897605531389218 0.86260350666634844 -0.24116554631735238 -0.81576638649106037
		 -2.0051156455114816 0.35748951030872234 -0.2684318973020936 -1.8450450571321859 -0.91797426133968474
		 0.26057726528001129 0.67536070616984833 -0.55769002176560245 -0.20931494635813694
		 -0.70231178477602396 -0.61679190927873084 1.2524426872206846 -0.09589870367790132 -0.30266797631558617
		 -0.1482835168096642 -0.68420599758732481 -1.4403841390245242 0.9755233777732315 -1.6540680007343191 -0.13357392534292034
		 0.36035343679910742 -0.46451859573662446 -1.0558131576663066 0.76167037849313601
		 -0.18529552145142603 -1.4026401511312878 0.0031560926963152269 0.57599056678971228 -1.3992855005317089
		 -1.0596729795767033 -0.98241926150592251 -1.721427931812024 0.5815940849193908 -0.37676682638106979 -0.65359161941981869
		 -1.0714614983654429 0.13474995764744557 0.072637639157216771 0.20555605620958772
		 -1.1969489913675637 -1.0343675034570701 -0.83856771514238637 -1.1282399989348313 -1.7596333196892822
		 -1.0183734176347572 -1.2020560335558437 -0.68715715056440485 0.062350634214078583
		 0.55726091014260182 -0.71926121881735905 1.2054013669242296 -0.29983522883834224 -1.4981987702733015
		 -0.11023274037720512 -0.39590994224241527 -0.66075696484186819 -0.72236652138128654
		 -0.02807875289644024 -0.6194702105226263 -1.5268366862619702 -0.46849478602384592 -0.29633813632066147
		 -0.23029280433811439 -2.0767778406211144 -1.9368136590634883 0.72043318764289099
		 -0.4812158010671489 -2.1456037538492918 0.36578227589421319 -1.1289137118987846 -1.187681792040717
		 -0.57843393648075225 -1.808991405210103 0.7669146079809146 0.081102383662005201 -0.43761173920137036 -0.49301553060301939
		 -0.75471295803828298 -0.46458526212459006 -0.09262355922146448 -0.27883268933491084
		 -0.81323005898155787 -0.14051870463277893 0.24179694659294457 -2.0798206912482078 -0.69046901005869399
		 1.397455452258231 -0.53566965179269699 0.46785373134818609 -0.85306470220014607 -0.87564399026055428 -1.8914025651353876
		 0.15750804522086811 0.083936797032722499 -0.4044719577910989 0.07518218848153993
		 0.5170699816328943 -1.5196850736479455 0.033393219599582837 0.56801906122785439 -0.19785676632411575
		 0.70771947867684737 -1.3407783593741081 0.29068778191745048 0.59298226442209623 -1.8029878207462708 -1.009362907908939
		 -0.20899375698693348 0.16363169628513319 -0.22276946049599855 -1.2500893844115126
		 -0.31677597051603712 -0.75323951136905332 -0.76008906177338531 -0.89333157270865304 0.48548276476361774
		 0.57275609448213349 -0.071574556455229921 -1.4359590843941981 -0.044376558108660903
		 0.16230854370654912 -0.63109875308423324 -0.54117134216552365 -0.57713547505083862 -0.89607489851810485
		 -0.47091880899250216 -2.1687602698351878 -0.93282175663529687 1.1475467411859597
		 -0.7851842531726938 -1.7257283542532769 -0.040292732982737792 -1.0278135073384569 -1.4909613763155045
		 -0.45490082460276543 -0.83754251241094058 -0.67234665689492701 1.2231295281409702
		 -1.23485012004434 -1.5006498540160913 -1.1771959153255012 -1.1321218302978628 -0.0507879295844752
		 1.0827160920665733 -1.163115323052248 -1.1216378534869016 -0.4437254035758656 -1.4122422752854256 -1.6122270541121049
		 0.12605910560408234 -1.1846356155584126 -0.36389855576871283 -0.055460238167023591
		 -0.11930732062948335 -1.5436717583109094 0.2974043025239041 -1.3573859976563805 -1.7138173548631521
		 0.46285172665067764 0.1996410595679059 -0.95919728171477847 -0.58709848182066371
		 -1.6393488688073186 ;
	setAttr ".id0" -type "doubleArray" 100 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
		 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43
		 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70
		 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97
		 98 99 ;
	setAttr ".nid" 100;
	setAttr ".nid0" 100;
	setAttr ".irbx" -type "string" "";
	setAttr ".irax" -type "string" "";
	setAttr ".icx" -type "string" "";
	setAttr ".cts" 1;
	setAttr ".prt" 7;
	setAttr ".boce" 0.89999997615814209;
	setAttr ".fron" 0.019999999552965164;
	setAttr ".cofl" 1;
	setAttr -s 2 ".fsc[0:1]"  0 1 1 1 0 1;
	setAttr -s 2 ".pfdo[0:1]"  0 1 1 1 0 1;
	setAttr ".vssc[0]"  0 1 1;
	setAttr ".stns[0]"  0 1 1;
	setAttr ".thr" 0;
	setAttr ".rdc[0]"  0 1 1;
	setAttr ".rci" 1;
	setAttr ".mssc[0]"  0 1 1;
	setAttr ".pfsc[0]"  0 1 1;
	setAttr ".frsc[0]"  0 1 1;
	setAttr ".stsc[0]"  0 1 1;
	setAttr ".clsc[0]"  0 1 1;
	setAttr ".bosc[0]"  0 1 1;
	setAttr ".opc[0]"  0 1 1;
	setAttr ".oci" 1;
	setAttr -s 2 ".cl";
	setAttr ".cl[0].clp" 0;
	setAttr ".cl[0].clc" -type "float3" 1 0 0 ;
	setAttr ".cl[0].cli" 1;
	setAttr ".cl[1].clp" 1;
	setAttr ".cl[1].clc" -type "float3" 0 0.5 1 ;
	setAttr ".cl[1].cli" 1;
	setAttr ".coi" 6;
	setAttr ".inca[0].incap" 0;
	setAttr ".inca[0].incac" -type "float3" 0 0 0 ;
	setAttr ".inca[0].incai" 1;
	setAttr -k on ".lifespan" 1;
createNode lightLinker -s -n "lightLinker1";
	rename -uid "F88686B9-4CAC-362B-EDC0-E19B5586E3A7";
	setAttr -s 3 ".lnk";
	setAttr -s 3 ".slnk";
createNode shapeEditorManager -n "shapeEditorManager";
	rename -uid "062F327D-411A-548F-D10E-ADB2D29ABC70";
createNode poseInterpolatorManager -n "poseInterpolatorManager";
	rename -uid "23F5E46F-4571-918A-7AD9-4EA2F7B6B4A0";
createNode displayLayerManager -n "layerManager";
	rename -uid "B8D4002C-43FE-926E-23E4-4197CECEB011";
createNode displayLayer -n "defaultLayer";
	rename -uid "5510B69C-43DF-F2C0-45E1-1DBEBE3271BD";
createNode renderLayerManager -n "renderLayerManager";
	rename -uid "3A1EFDFB-4E3B-2AFE-04B4-9088655DB899";
createNode renderLayer -n "defaultRenderLayer";
	rename -uid "042F3848-4B84-623C-FE19-CBA956B79E3C";
	setAttr ".g" yes;
createNode polyCube -n "polyCube1";
	rename -uid "5CF99807-4B3F-5C5B-5B74-1FA206D394C0";
	setAttr ".cuv" 4;
createNode deleteComponent -n "deleteComponent1";
	rename -uid "4E3B2146-495A-9555-D28D-0481ACD23E75";
	setAttr ".dc" -type "componentList" 1 "f[1]";
createNode shadingEngine -n "nParticleBallsSE";
	rename -uid "F650B898-4211-0CF8-DC7F-528F544E996C";
	setAttr ".ihi" 0;
	setAttr ".ro" yes;
createNode materialInfo -n "materialInfo1";
	rename -uid "D72CB7E9-4D54-C89D-FD49-23A876F13D58";
createNode particleSamplerInfo -n "particleSamplerInfo1";
	rename -uid "52BD14FB-4D23-11DA-7721-FB9BFA66F320";
createNode blinn -n "npBallsBlinn";
	rename -uid "48D5FC62-4D19-E7DA-47D0-2FA92ABE3EB9";
createNode particleCloud -n "npBallsVolume";
	rename -uid "36110146-4B44-6A70-A26D-159F0D2E4E11";
select -ne :time1;
	setAttr ".o" 1;
	setAttr ".unw" 1;
select -ne :hardwareRenderingGlobals;
	setAttr ".otfna" -type "stringArray" 22 "NURBS Curves" "NURBS Surfaces" "Polygons" "Subdiv Surface" "Particles" "Particle Instance" "Fluids" "Strokes" "Image Planes" "UI" "Lights" "Cameras" "Locators" "Joints" "IK Handles" "Deformers" "Motion Trails" "Components" "Hair Systems" "Follicles" "Misc. UI" "Ornaments"  ;
	setAttr ".otfva" -type "Int32Array" 22 0 1 1 1 1 1
		 1 1 1 0 0 0 0 0 0 0 0 0
		 0 0 0 0 ;
	setAttr ".fprt" yes;
select -ne :renderPartition;
	setAttr -s 3 ".st";
select -ne :renderGlobalsList1;
select -ne :defaultShaderList1;
	setAttr -s 6 ".s";
select -ne :postProcessList1;
	setAttr -s 2 ".p";
select -ne :defaultRenderingList1;
select -ne :initialShadingGroup;
	setAttr ".ro" yes;
select -ne :initialParticleSE;
	setAttr ".ro" yes;
select -ne :defaultRenderGlobals;
	setAttr ".ren" -type "string" "arnold";
select -ne :defaultResolution;
	setAttr ".pa" 1;
select -ne :hardwareRenderGlobals;
	setAttr ".ctrs" 256;
	setAttr ".btrs" 512;
select -ne :ikSystem;
	setAttr -s 4 ".sol";
connectAttr "deleteComponent1.og" "pCubeShape1.i";
connectAttr ":time1.o" "nucleus1.cti";
connectAttr "nParticleShape1.cust" "nucleus1.niao[0]";
connectAttr "nParticleShape1.stst" "nucleus1.nias[0]";
connectAttr ":time1.o" "nParticleShape1.cti";
connectAttr "nucleus1.noao[0]" "nParticleShape1.nxst";
connectAttr "nucleus1.stf" "nParticleShape1.stf";
connectAttr "nParticleShape1.incr" "nParticleShape1.rgbPP";
connectAttr "nParticleShape1.inor" "nParticleShape1.opacityPP";
connectAttr "nParticleShape1.inrr" "nParticleShape1.radiusPP";
relationship "link" ":lightLinker1" ":initialShadingGroup.message" ":defaultLightSet.message";
relationship "link" ":lightLinker1" ":initialParticleSE.message" ":defaultLightSet.message";
relationship "link" ":lightLinker1" "nParticleBallsSE.message" ":defaultLightSet.message";
relationship "shadowLink" ":lightLinker1" ":initialShadingGroup.message" ":defaultLightSet.message";
relationship "shadowLink" ":lightLinker1" ":initialParticleSE.message" ":defaultLightSet.message";
relationship "shadowLink" ":lightLinker1" "nParticleBallsSE.message" ":defaultLightSet.message";
connectAttr "layerManager.dli[0]" "defaultLayer.id";
connectAttr "renderLayerManager.rlmi[0]" "defaultRenderLayer.rlid";
connectAttr "polyCube1.out" "deleteComponent1.ig";
connectAttr "npBallsBlinn.oc" "nParticleBallsSE.ss";
connectAttr "npBallsVolume.oi" "nParticleBallsSE.vs";
connectAttr "nParticleShape1.iog" "nParticleBallsSE.dsm" -na;
connectAttr "nParticleBallsSE.msg" "materialInfo1.sg";
connectAttr "npBallsBlinn.msg" "materialInfo1.m";
connectAttr "particleSamplerInfo1.msg" "materialInfo1.t" -na;
connectAttr "particleSamplerInfo1.oc" "npBallsBlinn.c";
connectAttr "particleSamplerInfo1.ot" "npBallsBlinn.it";
connectAttr "particleSamplerInfo1.oi" "npBallsBlinn.ic";
connectAttr "particleSamplerInfo1.ot" "npBallsVolume.t";
connectAttr "particleSamplerInfo1.oc" "npBallsVolume.c";
connectAttr "particleSamplerInfo1.oi" "npBallsVolume.i";
connectAttr "nParticleBallsSE.pa" ":renderPartition.st" -na;
connectAttr "npBallsBlinn.msg" ":defaultShaderList1.s" -na;
connectAttr "npBallsVolume.msg" ":defaultShaderList1.s" -na;
connectAttr "defaultRenderLayer.msg" ":defaultRenderingList1.r" -na;
connectAttr "pCubeShape1.iog" ":initialShadingGroup.dsm" -na;
// End of ParticlesInCube.ma
