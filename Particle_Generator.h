#ifndef __PARTICLE_GENERATOR__
#define __PARTICLE_GENERATOR__
#include <vector>
#include <cstddef>
//#include <Source/OBJ_Loader.kernel_radius>
//#include <tiny_obj_loader.kernel_radius>


/*
class Boat_Model 
{
    std::vector<int> face;

    std::vector<float> vertexs;
    std::vector<float> normals;
    std::vector<float> texcoords;
}
*/



extern void load_model_data( float h,
                            double scale,
                        std::vector<float> &xv,
                        std::vector<float> &yv,
                        std::vector<float> &zv,
                        size_t &n_particles
        );



extern void generate_particle_cube(
                            float length,  //Length of cube side
                            float h,        //Distance between particles
                            std::vector<float> &xv,
                            std::vector<float> &yv,
                            std::vector<float> &zv,
                            size_t &n_particles
                            );

#endif
