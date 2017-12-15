//Maya ASCII 2017 scene
//Name: frame301.ma
//Last modified: Fri, Nov 24, 2017 05:44:03 PM
//Codeset 1252
requires maya "2017";currentUnit -l centimeter -a degree -t film;fileInfo "application" "maya";fileInfo "product" "Maya 2017";fileInfo "version" "2017";fileInfo "cutIdentifier" "201606150345-997974";fileInfo "osv" "Microsoft Windows 7 Business Edition, 64-bit Windows 7 Service Pack 1 (Build 7601)\n";fileInfo "license" "student";createNode transform -s -n "persp";	rename -uid "7185C08C-45B0-9914-5E7F-F58B088B941F";	setAttr ".v" no;	setAttr ".t" -type "double3" 28 21 28 ;	setAttr ".r" -type "double3" -27.938352729602379 44.999999999999972 -5.172681101354183e-014 ;createNode camera -s -n "perspShape" -p "persp";	rename -uid "B7B76270-45C3-73B5-C146-C9945484BF41";	setAttr -k off ".v" no;	setAttr ".fl" 34.999999999999993;	setAttr ".coi" 44.82186966202994;	setAttr ".imn" -type "string" "persp";	setAttr ".den" -type "string" "persp_depth";	setAttr ".man" -type "string" "persp_mask";	setAttr ".hc" -type "string" "viewSet -p %camera";	setAttr ".ai_translator" -type "string" "perspective";createNode transform -s -n "top";	rename -uid "AD93D0F2-479F-6FCC-C4A5-0582C3675BBD";	setAttr ".v" no;	setAttr ".t" -type "double3" 0 1000.1 0 ;	setAttr ".r" -type "double3" -89.999999999999986 0 0 ;createNode camera -s -n "topShape" -p "top";	rename -uid "AAFABFA7-4695-1258-6C75-8F93F3B1289B";	setAttr -k off ".v" no;	setAttr ".rnd" no;	setAttr ".coi" 1000.1;	setAttr ".ow" 30;	setAttr ".imn" -type "string" "top";	setAttr ".den" -type "string" "top_depth";	setAttr ".man" -type "string" "top_mask";	setAttr ".hc" -type "string" "viewSet -t %camera";	setAttr ".o" yes;	setAttr ".ai_translator" -type "string" "orthographic";createNode transform -s -n "front";	rename -uid "07423675-4AD0-D930-EFB8-D08CCD033B04";	setAttr ".v" no;	setAttr ".t" -type "double3" 0 0 1000.1 ;createNode camera -s -n "frontShape" -p "front";	rename -uid "9A7FC0F4-4EAD-EFA5-AD34-52A14EB2E365";	setAttr -k off ".v" no;	setAttr ".rnd" no;	setAttr ".coi" 1000.1;	setAttr ".ow" 30;	setAttr ".imn" -type "string" "front";	setAttr ".den" -type "string" "front_depth";	setAttr ".man" -type "string" "front_mask";	setAttr ".hc" -type "string" "viewSet -f %camera";	setAttr ".o" yes;	setAttr ".ai_translator" -type "string" "orthographic";createNode transform -s -n "side";	rename -uid "505A48E9-4E7B-864B-F0E8-BDA4AB8E4EE8";	setAttr ".v" no;	setAttr ".t" -type "double3" 1000.1 0 0 ;	setAttr ".r" -type "double3" 0 89.999999999999986 0 ;createNode camera -s -n "sideShape" -p "side";	rename -uid "D11E21BD-4D87-E472-ECBC-819FC310152D";	setAttr -k off ".v" no;	setAttr ".rnd" no;	setAttr ".coi" 1000.1;	setAttr ".ow" 30;	setAttr ".imn" -type "string" "side";	setAttr ".den" -type "string" "side_depth";	setAttr ".man" -type "string" "side_mask";	setAttr ".hc" -type "string" "viewSet -s %camera";	setAttr ".o" yes;	setAttr ".ai_translator" -type "string" "orthographic";createNode transform -n "pCube1";
	rename -uid "649B1D52-4142-B461-C245-F184D43846BB";
	setAttr ".s" -type "double3" 5 5 5;
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
	setAttr ".s" -type "double3" 5 5 5;
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
	setAttr ".pos0" -type "vectorArray" 69 0.0353294 0.0373149 0.0183657 6.29734e-05 0.0666168 0.0881607 0.0001 0.0433106 0.0489657 0.0867344 0.0639323 0.0277309 0.070343 0.0419991 0.0628279 0.0331629 0.0226382 0.079561 0.0411952 0.0860848 0.122382 0.131182 0.0350643 0.0938578 8.45392e-05 0.0942518 0.151726 0.0697309 0.00101431 8.77055e-06 0.0279572 0.0299286 0.16091 0.0879154 0.0793201 0.134819 0.0864495 2.54579e-05 0.0763604 0.0950569 2.00872e-05 0.121552 0.171337 0.0598685 0.00645337 0.187265 0.0533355 0.075441 0.189648 0.00690748 0.100713 0.00452473 0.0508312 0.182788 0.0564146 0.0368194 0.217975 0.168339 0.013892 0.166791 0.00342444 0.0833428 0.056042 0.0596094 0.0662751 0.187692 0.105531 0.0278866 0.184126 0.114615 0.0647397 0.184345 0.0528134 8.15703e-05 0.149642 0.138073 0.0971859 0.0639747 0.166425 3.48078e-05 7.43397e-05 0.216873 7.7389e-05 1.01742e-05 0.140493 8.50195e-05 0.077296 0.169625 0.055521 0.197917 0.135215 0.043275 0.0497492 0.128155 0.114878 0.138626 0.0755599 0.0335348 0.117373 0.0404044 0.0001 0.210353 0.0503541 0.098483 0.0148898 0.124805 0.0001 0.218188 0.150873 0.0296309 0.225672 0.121243 5.73758e-05 8.4238e-06 0.172744 6.52277e-05 0.04708 0.196442 0.0295849 0.0201536 0.225297 0.0306996 0.0759444 8.97402e-05 1.13821e-05 0.224349 0.20018 0.0552912 0.129298 0.109752 0.0142055 0.0414867 0.054385 0.0001 0.0482098 0.229869 1.42283e-05 0.0504225 0.141684 0.0465591 9.33315e-05 7.83409e-05 6.42838e-06 0.0987615 0.107373 0.0967 0.204509 0.201505 1.15495e-05 0.149564 0.0789439 8.89005e-05 0.172715 0.182641 3.72873e-05 0.221563 0.0514369 3.90958e-05 0.107364 6.40019e-05 0.0361653 0.122519 0.247378 0.0377447 0.127364 0.0236992 0.0001 0.0001 0.0843836 0.0749733 0.0802116 0.231257 0.0001 0.187132 0.0850734 0.0408106 6.27914e-05 0.0853663 9.71026e-05 0.242428 9.51537e-05 0.0001 0.0366522 0.141097 9.21722e-05 0.127371 0.211511 0.0318375 0.190036 0.148151 0.0447878 0.135803 0.200174 0.0989045 0.100078 5.50446e-05 0.0354644 8.96233e-05 9.08735e-05 6.37756e-05 0.154465 0.207282 0.0907983 0.0427857 0.221238 0.102288 0.164118 ;
	setAttr ".id0" -type "doubleArray" 69 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 ;
	setAttr ".nid" 69;
	setAttr ".nid0" 69;
	setAttr ".irbx" -type "string" "";	setAttr ".irax" -type "string" "";	setAttr ".icx" -type "string" "";	setAttr ".cts" 1;	setAttr ".prt" 7;	setAttr ".boce" 0.89999997615814209;	setAttr ".fron" 0.019999999552965164;	setAttr ".cofl" 1;	setAttr -s 2 ".fsc[0:1]"  0 1 1 1 0 1;	setAttr -s 2 ".pfdo[0:1]"  0 1 1 1 0 1;	setAttr ".vssc[0]"  0 1 1;	setAttr ".stns[0]"  0 1 1;	setAttr ".thr" 0;	setAttr ".rdc[0]"  0 1 1;	setAttr ".rci" 1;	setAttr ".mssc[0]"  0 1 1;	setAttr ".pfsc[0]"  0 1 1;	setAttr ".frsc[0]"  0 1 1;	setAttr ".stsc[0]"  0 1 1;	setAttr ".clsc[0]"  0 1 1;	setAttr ".bosc[0]"  0 1 1;	setAttr ".opc[0]"  0 1 1;	setAttr ".oci" 1;	setAttr -s 2 ".cl";	setAttr ".cl[0].clp" 0;	setAttr ".cl[0].clc" -type "float3" 1 0 0 ;	setAttr ".cl[0].cli" 1;	setAttr ".cl[1].clp" 1;	setAttr ".cl[1].clc" -type "float3" 0 0.5 1 ;	setAttr ".cl[1].cli" 1;	setAttr ".coi" 6;	setAttr ".inca[0].incap" 0;	setAttr ".inca[0].incac" -type "float3" 0 0 0 ;	setAttr ".inca[0].incai" 1;	setAttr -k on ".lifespan" 1;createNode lightLinker -s -n "lightLinker1";	rename -uid "F88686B9-4CAC-362B-EDC0-E19B5586E3A7";	setAttr -s 3 ".lnk";	setAttr -s 3 ".slnk";createNode shapeEditorManager -n "shapeEditorManager";	rename -uid "062F327D-411A-548F-D10E-ADB2D29ABC70";createNode poseInterpolatorManager -n "poseInterpolatorManager";	rename -uid "23F5E46F-4571-918A-7AD9-4EA2F7B6B4A0";createNode displayLayerManager -n "layerManager";	rename -uid "B8D4002C-43FE-926E-23E4-4197CECEB011";createNode displayLayer -n "defaultLayer";	rename -uid "5510B69C-43DF-F2C0-45E1-1DBEBE3271BD";createNode renderLayerManager -n "renderLayerManager";	rename -uid "3A1EFDFB-4E3B-2AFE-04B4-9088655DB899";createNode renderLayer -n "defaultRenderLayer";	rename -uid "042F3848-4B84-623C-FE19-CBA956B79E3C";	setAttr ".g" yes;createNode polyCube -n "polyCube1";	rename -uid "5CF99807-4B3F-5C5B-5B74-1FA206D394C0";	setAttr ".cuv" 4;createNode deleteComponent -n "deleteComponent1";	rename -uid "4E3B2146-495A-9555-D28D-0481ACD23E75";	setAttr ".dc" -type "componentList" 1 "f[1]";createNode shadingEngine -n "nParticleBallsSE";	rename -uid "F650B898-4211-0CF8-DC7F-528F544E996C";	setAttr ".ihi" 0;	setAttr ".ro" yes;createNode materialInfo -n "materialInfo1";	rename -uid "D72CB7E9-4D54-C89D-FD49-23A876F13D58";createNode particleSamplerInfo -n "particleSamplerInfo1";	rename -uid "52BD14FB-4D23-11DA-7721-FB9BFA66F320";createNode blinn -n "npBallsBlinn";	rename -uid "48D5FC62-4D19-E7DA-47D0-2FA92ABE3EB9";createNode particleCloud -n "npBallsVolume";	rename -uid "36110146-4B44-6A70-A26D-159F0D2E4E11";select -ne :time1;	setAttr ".o" 1;	setAttr ".unw" 1;select -ne :hardwareRenderingGlobals;	setAttr ".otfna" -type "stringArray" 22 "NURBS Curves" "NURBS Surfaces" "Polygons" "Subdiv Surface" "Particles" "Particle Instance" "Fluids" "Strokes" "Image Planes" "UI" "Lights" "Cameras" "Locators" "Joints" "IK Handles" "Deformers" "Motion Trails" "Components" "Hair Systems" "Follicles" "Misc. UI" "Ornaments"  ;	setAttr ".otfva" -type "Int32Array" 22 0 1 1 1 1 1		 1 1 1 0 0 0 0 0 0 0 0 0		 0 0 0 0 ;	setAttr ".fprt" yes;select -ne :renderPartition;	setAttr -s 3 ".st";select -ne :renderGlobalsList1;select -ne :defaultShaderList1;	setAttr -s 6 ".s";select -ne :postProcessList1;	setAttr -s 2 ".p";select -ne :defaultRenderingList1;select -ne :initialShadingGroup;	setAttr ".ro" yes;select -ne :initialParticleSE;	setAttr ".ro" yes;select -ne :defaultRenderGlobals;	setAttr ".ren" -type "string" "arnold";select -ne :defaultResolution;	setAttr ".pa" 1;select -ne :hardwareRenderGlobals;	setAttr ".ctrs" 256;	setAttr ".btrs" 512;select -ne :ikSystem;	setAttr -s 4 ".sol";connectAttr "deleteComponent1.og" "pCubeShape1.i";connectAttr ":time1.o" "nucleus1.cti";connectAttr "nParticleShape1.cust" "nucleus1.niao[0]";connectAttr "nParticleShape1.stst" "nucleus1.nias[0]";connectAttr ":time1.o" "nParticleShape1.cti";connectAttr "nucleus1.noao[0]" "nParticleShape1.nxst";connectAttr "nucleus1.stf" "nParticleShape1.stf";connectAttr "nParticleShape1.incr" "nParticleShape1.rgbPP";connectAttr "nParticleShape1.inor" "nParticleShape1.opacityPP";connectAttr "nParticleShape1.inrr" "nParticleShape1.radiusPP";relationship "link" ":lightLinker1" ":initialShadingGroup.message" ":defaultLightSet.message";relationship "link" ":lightLinker1" ":initialParticleSE.message" ":defaultLightSet.message";relationship "link" ":lightLinker1" "nParticleBallsSE.message" ":defaultLightSet.message";relationship "shadowLink" ":lightLinker1" ":initialShadingGroup.message" ":defaultLightSet.message";relationship "shadowLink" ":lightLinker1" ":initialParticleSE.message" ":defaultLightSet.message";relationship "shadowLink" ":lightLinker1" "nParticleBallsSE.message" ":defaultLightSet.message";connectAttr "layerManager.dli[0]" "defaultLayer.id";connectAttr "renderLayerManager.rlmi[0]" "defaultRenderLayer.rlid";connectAttr "polyCube1.out" "deleteComponent1.ig";connectAttr "npBallsBlinn.oc" "nParticleBallsSE.ss";connectAttr "npBallsVolume.oi" "nParticleBallsSE.vs";connectAttr "nParticleShape1.iog" "nParticleBallsSE.dsm" -na;connectAttr "nParticleBallsSE.msg" "materialInfo1.sg";connectAttr "npBallsBlinn.msg" "materialInfo1.m";connectAttr "particleSamplerInfo1.msg" "materialInfo1.t" -na;connectAttr "particleSamplerInfo1.oc" "npBallsBlinn.c";connectAttr "particleSamplerInfo1.ot" "npBallsBlinn.it";connectAttr "particleSamplerInfo1.oi" "npBallsBlinn.ic";connectAttr "particleSamplerInfo1.ot" "npBallsVolume.t";connectAttr "particleSamplerInfo1.oc" "npBallsVolume.c";connectAttr "particleSamplerInfo1.oi" "npBallsVolume.i";connectAttr "nParticleBallsSE.pa" ":renderPartition.st" -na;connectAttr "npBallsBlinn.msg" ":defaultShaderList1.s" -na;connectAttr "npBallsVolume.msg" ":defaultShaderList1.s" -na;connectAttr "defaultRenderLayer.msg" ":defaultRenderingList1.r" -na;connectAttr "pCubeShape1.iog" ":initialShadingGroup.dsm" -na;
// End of frame301.ma;