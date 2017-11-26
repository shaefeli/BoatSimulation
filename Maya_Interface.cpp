#include <iostream>
#include <fstream>
#include <sstream>
#include "Maya_Interface.h"
#include <string.h>
using namespace std;

Maya_Interface::Maya_Interface(size_t nr_particles){
	//fileStart = read file with start
	ifstream t1("maya_utilities/typicalFileStartMaya");
	ifstream t2("maya_utilities/typicalFileEndMaya");
	string temp = "";
	string temp2 = "";
	while(getline(t1,temp)){
		temp2 = temp2+temp;
	}
	Maya_Interface::fileStart = temp2;
	temp="";
	temp2="";
	while(getline(t2,temp)){
		temp2=temp2+temp;
	}
	Maya_Interface::fileEnd = temp2;
	
	//initialize "characteristics that are not position but still have to go to maya file" 
}
void Maya_Interface::writeToMaya(size_t frameNr, float* interleaved_positions, size_t nr_particles){
	ofstream mayaFile;
	stringstream header;
	header << "frame"<<frameNr<<".ma";
  	mayaFile.open(header.str().c_str());
  	
  	//Header
  	header.str(string());
  	header <<"//Maya ASCII 2017 scene\n" <<  "Name: frame" << frameNr << ".ma\n"<<"//Last modified: Fri, Nov 24, 2017 05:44:03 PM\n"<<"//Codeset 1252\n";
    double cube_sizeX = 10;
    double cube_sizeY = 2;
    double cube_sizeZ = 10;
    stringstream cube;
  	string cubeInfos1 = string("createNode transform -n \"pCube1\";\n")+
							string("rename -uid \"649B1D52-4142-B461-C245-F184D43846BB\";\n");
	string cubeInfos2 = string("createNode mesh -n \"pCubeShape1\" -p \"pCube1\";\n")+			
							string("rename -uid \"60008B7B-4CAB-75B2-C86E-43A99298F542\";\n")+
							string("setAttr -k off \".v\";\n")+
							string("setAttr \".vir\" yes;\n")+
							string("setAttr \".vif\" yes;\n")+
							string("setAttr \".pv\" -type \"double2\" 0.5 0.375 ;\n")+
							string("setAttr \".uvst[0].uvsn\" -type \"string\" \"map1\";\n")+
							string("setAttr \".cuvs\" -type \"string\" \"map1\";\n")+
							string("setAttr \".dcc\" -type \"string\" \"Ambient+Diffuse\";\n")+
							string("setAttr \".covm[0]\"  0 1 1;\n")+
							string("setAttr \".cdvm[0]\"  0 1 1;\n")+
							string("setAttr \".ai_translator\" -type \"string\" \"polymesh\";\n");
	cube << cubeInfos1 << "setAttr \".s\" -type \"double3\" " << cube_sizeX << " " << cube_sizeY << " " << cube_sizeZ << ";\n" << cubeInfos2;
	//Question: how to access fields from constructor? 
	//cout << header.str() << Maya_Interface::fileStart << cube.str() << Maya_Interface::fileEnd << "// End of frame" << frameNr << ".ma" << endl;
  	string particles = "";

  	//Write to file
  	mayaFile << "Buongiorno";
  	mayaFile.close();
}