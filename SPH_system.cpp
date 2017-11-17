#include "SPH_system.h"

SPH_system::SPH_system(){
	//Construct the initial state of the particles (this may be done later with an obj file)
	Vector3T<double> pos1 = Vector3T<double>(0.5,0.5,1);
	Vector3T<double> pos2 = Vector3T<double>(0.5,0.5,1);
	position[0] = pos1;
	position[1] = pos2;
	signs[0]=1;
	signs[1]=1;
}

//Main function that updates our particle system
void SPH_system::update(){
	double posX1 = this->position[0][0];
	double posY2 = this->position[1][1];

	if(posX1<=0 || posX1>=2.0){
		signs[0] = signs[0]*-1;
	}
	if(posY2<=0 || posY2>=2.0){
		signs[1] = signs[1]*-1;
	}
	this->position[0][0] += signs[0]*0.001;
	this->position[1][1] += signs[1]*0.001;
}
