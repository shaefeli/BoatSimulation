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
	
	stringstream idis;
	idis << "\tsetAttr \".id0\" -type \"doubleArray\" " << nr_particles << " ";
	for(int i=0; i<nr_particles;i++){
		idis << i << " ";
	}
	idis << ";\n";
	Maya_Interface::ids = idis.str();
}
void Maya_Interface::writeToMaya(size_t frameNr, float* x ,float* y, float* z, size_t nr_particles){

	//Open the file to write
	ofstream mayaFile;
	stringstream header;
	header << "frameFiles/frame"<<frameNr<<".ma";
  	mayaFile.open(header.str().c_str());
  	
  	
  	//Header
  	header.str(string());
  	header <<"//Maya ASCII 2017 scene\n" <<  "//Name: frame" << frameNr << ".ma\n"<<"//Last modified: Fri, Nov 24, 2017 05:44:03 PM\n"<<"//Codeset 1252\n";

  	//Cube infos
    double cube_sizeX = 5;
    double cube_sizeY = 5;
    double cube_sizeZ = 5;
    stringstream cube;
  	string cubeInfos1 = string("createNode transform -n \"pCube1\";\n")+
							string("\trename -uid \"649B1D52-4142-B461-C245-F184D43846BB\";\n");
	string cubeInfos2 = string("createNode mesh -n \"pCubeShape1\" -p \"pCube1\";\n")+			
							string("\trename -uid \"60008B7B-4CAB-75B2-C86E-43A99298F542\";\n")+
							string("\tsetAttr -k off \".v\";\n")+
							string("\tsetAttr \".vir\" yes;\n")+
							string("\tsetAttr \".vif\" yes;\n")+
							string("\tsetAttr \".pv\" -type \"double2\" 0.5 0.375 ;\n")+
							string("\tsetAttr \".uvst[0].uvsn\" -type \"string\" \"map1\";\n")+
							string("\tsetAttr \".cuvs\" -type \"string\" \"map1\";\n")+
							string("\tsetAttr \".dcc\" -type \"string\" \"Ambient+Diffuse\";\n")+
							string("\tsetAttr \".covm[0]\"  0 1 1;\n")+
							string("\tsetAttr \".cdvm[0]\"  0 1 1;\n")+
							string("\tsetAttr \".ai_translator\" -type \"string\" \"polymesh\";\n");
	cube << cubeInfos1 << "\tsetAttr \".s\" -type \"double3\" " << cube_sizeX << " " << cube_sizeY << " " << cube_sizeZ << ";\n" << cubeInfos2;


  	//Particle infos
  	double partScaleX = 5;
  	double partScaleY = 5;
  	double partScaleZ = 5;
  	stringstream partScale;
  	partScale << "\tsetAttr \".s\" -type \"double3\" " << partScaleX << " " << partScaleY << " " << partScaleZ << ";\n";
  	string partInfo1 = string("createNode transform -n \"nParticle1\";\n")+
						string("\trename -uid \"287A3003-4FDB-C9D2-628D-23864FCF4CEC\";\n");
	string partInfo2 =	string("createNode nParticle -n \"nParticleShape1\" -p \"nParticle1\";\n")+
							string("\trename -uid \"B44B925C-457C-8339-4764-FBA63FC2DCFD\";\n")+
							string("\taddAttr -s false -ci true -sn \"lifespanPP\" -ln \"lifespanPP\" -dt \"doubleArray\";\n")+
							string("\taddAttr -ci true -h true -sn \"lifespanPP0\" -ln \"lifespanPP0\" -dt \"doubleArray\";\n")+
							string("\taddAttr -ci true -sn \"lifespan\" -ln \"lifespan\" -at \"double\";\n")+
							string("\taddAttr -s false -ci true -sn \"rgbPP\" -ln \"rgbPP\" -dt \"vectorArray\";\n")+
							string("\taddAttr -ci true -h true -sn \"rgbPP0\" -ln \"rgbPP0\" -dt \"vectorArray\";\n")+
							string("\taddAttr -s false -ci true -sn \"opacityPP\" -ln \"opacityPP\" -dt \"doubleArray\";\n")+
							string("\taddAttr -ci true -h true -sn \"opacityPP0\" -ln \"opacityPP0\" -dt \"doubleArray\";\n")+
							string("\taddAttr -s false -ci true -sn \"radiusPP\" -ln \"radiusPP\" -dt \"doubleArray\";\n")+
							string("\taddAttr -ci true -h true -sn \"radiusPP0\" -ln \"radiusPP0\" -dt \"doubleArray\";\n")+
							string("\tsetAttr -k off \".v\";\n")+
							string("\tsetAttr \".gf\" -type \"Int32Array\" 0 ;\n");
	stringstream particuleInfos;
	particuleInfos << "createNode nucleus -n \"nucleus1\";\n" << "\trename -uid \"C11D1C29-44B6-A668-506C-24B574D4BE7C\";\n" << partInfo1 << partScale.str() << partInfo2;


	//particle positions
  	stringstream positions;
  	positions << "\tsetAttr \".pos0\" -type \"vectorArray\" " << nr_particles << " ";
  	for (int i=0; i<nr_particles; i++){
  		positions  << x[i] << " " << y[i] << " " << z[i] << " ";
  	}
  	positions << ";\n";


  	//Before end
  	stringstream afterParticles;
  	afterParticles << "\tsetAttr \".nid\" " << nr_particles <<";\n" << "\tsetAttr \".nid0\" " << nr_particles <<";\n";

  	//Write everything to file
  	mayaFile << header.str() << this->fileStart << cube.str() << particuleInfos.str() << positions.str() << this->ids << afterParticles.str() << this->fileEnd << "\n// End of frame" << frameNr << ".ma;"; 
  	mayaFile.close();
}