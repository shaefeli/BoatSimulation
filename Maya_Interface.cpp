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
	
	stringstream ids;
	ids << "setAttr \".id0\" -type \"doubleArray\" " << nr_particles << " ";
	for(int i=0; i<nr_particles;i++){
		ids << i << " ";
	}
	ids << ";\n";
}
void Maya_Interface::writeToMaya(size_t frameNr, float* interleaved_positions, size_t nr_particles){

	//Open the file to write
	ofstream mayaFile;
	stringstream header;
	header << "frameFiles/frame"<<frameNr<<".ma";
  	mayaFile.open(header.str().c_str());
  	
  	
  	//Header
  	header.str(string());
  	header <<"//Maya ASCII 2017 scene\n" <<  "//Name: frame" << frameNr << ".ma\n"<<"//Last modified: Fri, Nov 24, 2017 05:44:03 PM\n"<<"//Codeset 1252\n";

  	//Cube infos
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


  	//Particle infos
  	double partScaleX = 1;
  	double partScaleY = 1;
  	double partScaleZ = 1;
  	string partInfo1 = string("createNode transform -n \"nParticle1\";\n")+
						string("rename -uid \"287A3003-4FDB-C9D2-628D-23864FCF4CEC\";\n");

	string partInfo2 =	string("createNode nParticle -n \"nParticleShape1\" -p \"nParticle1\";\n")+
							string("rename -uid \"B44B925C-457C-8339-4764-FBA63FC2DCFD\";\n")+
							string("addAttr -s false -ci true -sn \"lifespanPP\" -ln \"lifespanPP\" -dt \"doubleArray\";\n")+
							string("addAttr -ci true -h true -sn \"lifespanPP0\" -ln \"lifespanPP0\" -dt \"doubleArray\";\n")+
							string("addAttr -ci true -sn \"lifespan\" -ln \"lifespan\" -at \"double\";\n")+
							string("addAttr -s false -ci true -sn \"rgbPP\" -ln \"rgbPP\" -dt \"vectorArray\";\n")+
							string("addAttr -ci true -h true -sn \"rgbPP0\" -ln \"rgbPP0\" -dt \"vectorArray\";\n")+
							string("addAttr -s false -ci true -sn \"opacityPP\" -ln \"opacityPP\" -dt \"doubleArray\";\n")+
							string("addAttr -ci true -h true -sn \"opacityPP0\" -ln \"opacityPP0\" -dt \"doubleArray\";\n")+
							string("addAttr -s false -ci true -sn \"radiusPP\" -ln \"radiusPP\" -dt \"doubleArray\";\n")+
							string("addAttr -ci true -h true -sn \"radiusPP0\" -ln \"radiusPP0\" -dt \"doubleArray\";\n")+
							string("setAttr -k off \".v\";\n")+
							string("setAttr \".gf\" -type \"Int32Array\" 0 ;\n");
	stringstream particuleInfos;
	particuleInfos << partInfo1 << "setAttr \".s\" -type \"double3\" " << partScaleX << " " << partScaleY << " " << partScaleZ << ";\n" << partInfo2;


	//particle positions
  	stringstream positions;
  	positions << "setAttr \".pos0\" -type \"vectorArray\" " << nr_particles << " ";
  	for (int i=0; i<nr_particles; i++){
  		positions  << interleaved_positions[i*3] << " " << interleaved_positions[i*3+1] << " " << interleaved_positions[i*3+2] << " ";
  	}
  	positions << ";\n";


  	//Before end
  	stringstream afterParticles;
  	afterParticles << "createNode nucleus -n \"nucleus1\";\n" << "rename -uid \"C11D1C29-44B6-A668-506C-24B574D4BE7C\";\n" << "setAttr \".nid\" " << nr_particles <<";\n" << "setAttr \".nid0\" " << nr_particles <<";\n";

  	//Write everything to file
  	mayaFile << header.str() << "fileStart\n" << cube.str() << particuleInfos.str() << positions.str() << "\nids\n" << afterParticles.str() << "\nfileEnd\n" << "// End of frame" << frameNr << ".ma;"; 
  	mayaFile.close();
}