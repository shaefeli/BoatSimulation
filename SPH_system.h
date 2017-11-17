#include "Vector3T.h"

#define nr_particles 2

class SPH_system{
public:
	//Here have to go every field that are of interest (position, velocity, pressure,...)
	//For the moment only the position and a dummy for displaying purpose
	Vector3T<double> position [nr_particles];
	int signs [nr_particles];
public:
	SPH_system();
	void update();
};