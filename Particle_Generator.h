#ifndef __PARTICLE_GENERATOR__
#define __PARTICLE_GENERATOR__
#include <vector>
#include <cstddef>
//#include <Source/OBJ_Loader.h>
//#include <tiny_obj_loader.h>


/*
class Boat_Model 
{
    std::vector<int> face;

    std::vector<float> vertexs;
    std::vector<float> normals;
    std::vector<float> texcoords;
}
*/



extern void load_model_data();



extern void generate_particle_cube(
                            float length,  //Length of cube side
                            float h,        //Distance between particles
                            std::vector<float> &xv,
                            std::vector<float> &yv,
                            std::vector<float> &zv,
                            size_t &n_particles 
                            );

#endif
