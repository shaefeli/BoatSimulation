#include <cstddef>
#include "Uniform_Grid.h"
#define nr_particles 2

#ifndef _PBS_BASIC_SPH_SYSTEM_
#define _PBS_BASIC_SPH_SYSTEM_

typedef struct {
    size_t n_particles;
    //Position
    float *x;
    float *y;
    float *z;

    //Velocity
    float *vx;
    float *vy;
    float *vz;

} particle_information_t;

class Basic_SPH_System{
    private:
        //Update all the particles
        void update_velocities_dummy(float dt);
        void update_positions_dummy(float dt);
        
        //Update one particle in particular
        void update_particle_position_dummy(int i, float dt);
        void update_particle_velocity_dummy(int i, float dt);

        size_t get_particle_number();

        float b_min_x, b_min_y, b_min_z;
        float b_max_x, b_max_y, b_max_z;
    public:
        //Neighbor-search-structure
        Uniform_Grid uniform_grid;

        //Container for the particle data
        particle_information_t particles;

        Basic_SPH_System(   size_t n_particles,
                            float b_min_x,          //Boundary values
                            float b_min_y,
                            float b_min_z,
                            float b_max_x,
                            float b_max_y,
                            float b_max_z,
                            float cell_x,
                            float cell_y,
                            float cell_z
                          );
 
        //run one step of the simulation
        void run_step(float dt);

};

#endif
