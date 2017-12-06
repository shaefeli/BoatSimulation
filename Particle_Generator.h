#ifndef __PARTICLE_GENERATOR__
#define __PARTICLE_GENERATOR__
#include <vector>
#include <cstddef>



extern void generate_particle_cube(
                            float length,  //Length of cube side
                            float h,        //Distance between particles
                            std::vector<float> &xv,
                            std::vector<float> &yv,
                            std::vector<float> &zv,
                            size_t &n_particles 
                            );

#endif
