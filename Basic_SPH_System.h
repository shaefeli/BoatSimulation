#include <cstddef>
#include "Uniform_Grid.h"
#include "Vector3T.h"

#define nr_particles 2

#ifndef _PBS_BASIC_SPH_SYSTEM_
#define _PBS_BASIC_SPH_SYSTEM_

typedef struct {
    float dt;
    float h;
    float rho0;
    float k;
    float mu;
    float g;
} SimState;


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
    
    // !Total force acting
    float *Fx;
    float *Fy;
    float *Fz;
    
    // densities
    float *rho;

    // pressures
    float *p;

    float mass; // for now all particles have the same mass
    
} particle_information_t;

class Basic_SPH_System{
private:
    SimState simState;
    
    //Update all the particles
    void update_velocities_dummy(float dt);
    void update_positions_dummy(float dt);

    //Update one particle in particular
    void update_particle_position_dummy(int i, float dt);
    void update_particle_velocity_dummy(int i, float dt);



    size_t get_particle_number();

    float b_min_x, b_min_y, b_min_z;
    float b_max_x, b_max_y, b_max_z;
    
    /**
        Calc fields and fill appropriate arrays
     */
    void calculate_Forces();
    void calculate_Pressures();
//    void calculate_Fvisc();
    void calculate_Densities();
    
    /**
     Different kernels for different properties
     */
    float evalKernel_poly6(int &i, int &j, float &h);
    float evalKernel_spiky(int i, int j, float h);
    float evalkernel_spiky_gradient(int &i, int &j, float &h);

    float evalKernel_visc( int i, int j, float h);
    float evalKernel_visc_laplacian(int &i, int &j, float &h);
    
    // distance between particles
    Vector3T<float> ij_vector(int i,int j);
    float distanceIJ(int i, int j);
    
    float h9,h6, h3, h2;
    float h9_315; // constants for kernels
    float h6_15;
    float h6_15_grad;
    float h3_15;
    float h3_15_visc;

    bool flag_finilized = false;
    
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

    void setSimState( SimState state ); // set sim parameters
    void finilizeInit(); // check if all parameters are provided and the set the flag to true
    
    //run one step of the simulation ( could adjust time step )
    void run_step(float dt);
    ~Basic_SPH_System();
};

#endif
