#include <iostream>
#include <fstream>
#include "Maya_Interface.h"
using namespace std;

void Maya_Interface::init(size_t nr_particles){
	//fileStart = read file with start
	ifstream t1("maya_utilities/typicalFileStartMaya")
	ifstream t2("maya_utilities/typicalFileEndMaya")

	//Create string that creates "useless" feaure values of particles
	fileStart << t1.rdbuf();
	fileEnd << t2.rdbuf();
	
	//initialize characteristics 
}
void Maya_Interface::writeToMaya(size_t frameNr, float* interleaved_positions, size_t nr_particles){
	ofstream mayaFile;
  	mayaFile.open ("frame"+frameNr+".ma");

  	//Header
  	string lastModified = "//Last modified: Fri, Nov 24, 2017 05:44:03 PM" 
  	string fileName = "Name: frame"+frameNr+".ma";
  	string Header = "//Maya ASCII 2017 scene\n"+
  						fileName+"\n"+
  						lastModified+"\n"+
  						"//Codeset 1252";

    double cube_sizeX = 10;
    //this is the height (y)
    double cube_sizeY = 2;
    double cube_sizeZ = 10;
  	string cubeInfos = "createNode transform -n \"pCube1\";\n
				rename -uid \"649B1D52-4142-B461-C245-F184D43846BB\";\n
				setAttr \".s\" -type \"double3\" "+cube_sizeX+" "+cube_sizeY+" "+ cube_sizeZ+";\n
				createNode mesh -n \"pCubeShape1\" -p \"pCube1\";\n
				rename -uid \"60008B7B-4CAB-75B2-C86E-43A99298F542\";\n
				setAttr -k off \".v\";\n
				setAttr \".vir\" yes;\n
				setAttr \".vif\" yes;\n
				setAttr \".pv\" -type \"double2\" 0.5 0.375 ;\n
				setAttr \".uvst[0].uvsn\" -type \"string\" \"map1\";\n
				setAttr \".cuvs\" -type \"string\" \"map1\";\n
				setAttr \".dcc\" -type \"string\" \"Ambient+Diffuse\";\n
				setAttr \".covm[0]\"  0 1 1;\n
				setAttr \".cdvm[0]\"  0 1 1;\n
				setAttr \".ai_translator\" -type \"string\" \"polymesh\";\n"


  	string particles = "";

  	string endLine = "// End of frame"+frameNr+".ma\n";


  	mayaFile << Header+fileStart.str()+cubeInfos+particles+fileEnd.str()+endLine;
  	mayaFile.close();
}