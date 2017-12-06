#include "Particle_Generator.h"


//Object centered around the origin
void generate_particle_cube(
                            float length,  //Length of cube side
                            float h,        //Distance between particles
                            std::vector<float> &xv,
                            std::vector<float> &yv,
                            std::vector<float> &zv,
                            size_t &n_particles 
                            )
{
    float xmin = -length;
    float ymin = -length;
    float zmin = -length;
    float xmax =  length;
    float ymax =  length;
    float zmax =  length;
    
    size_t x_particles = (xmax-xmin)/h;
    size_t y_particles = (ymax-ymin)/h;
    size_t z_particles = (zmax-zmin)/h;
    xv = std::vector<float>(x_particles*y_particles*z_particles);
    yv = std::vector<float>(x_particles*y_particles*z_particles);
    zv = std::vector<float>(x_particles*y_particles*z_particles);
    for( size_t i = 0; i < x_particles; i++ ) {
        for( size_t j = 0; j < y_particles; j++ ) {
            for( size_t k = 0; k < z_particles; k++ ) {
                xv[i*y_particles*z_particles + j*z_particles + k] = i*h+xmin;
                yv[i*y_particles*z_particles + j*z_particles + k] = j*h+ymin;
                zv[i*y_particles*z_particles + j*z_particles + k] = k*h+zmin;
            }
        }
    }
    n_particles = x_particles*y_particles*z_particles;
}

