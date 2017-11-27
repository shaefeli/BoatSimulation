#include "Basic_SPH_System.h"
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <math.h>

Basic_SPH_System::Basic_SPH_System( size_t n_particles,
                                    float b_min_x,//Boundary values
                                    float b_min_y,
                                    float b_min_z,
                                    float b_max_x,
                                    float b_max_y,
                                    float b_max_z,
                                    float cell_x,
                                    float cell_y,
                                    float cell_z
                                  ) :
    uniform_grid(b_min_x,b_min_y,b_min_z,
                 b_max_x,b_max_y,b_max_z,
                 cell_x,cell_y,cell_z),
    b_min_x(b_min_x),b_min_y(b_min_y),b_min_z(b_min_z),
    b_max_x(b_max_x),b_max_y(b_max_y),b_max_z(b_max_z)
{
    particles.n_particles = n_particles;

    particles.x = (float *)malloc(n_particles*sizeof(float));
    particles.y = (float *)malloc(n_particles*sizeof(float));
    particles.z = (float *)malloc(n_particles*sizeof(float));

    particles.vx = (float *)malloc(n_particles*sizeof(float));
    particles.vy = (float *)malloc(n_particles*sizeof(float));
    particles.vz = (float *)malloc(n_particles*sizeof(float));
    
    particles.Fx = (float *)malloc(n_particles*sizeof(float));
    particles.Fy = (float *)malloc(n_particles*sizeof(float));
    particles.Fz = (float *)malloc(n_particles*sizeof(float));
    
    particles.rho = (float *)malloc(n_particles*sizeof(float));

    for( size_t i = 0; i < n_particles; i++ ) {
       particles.x[i] = (float(rand())/RAND_MAX);
       particles.y[i] = (float(rand())/RAND_MAX);
       particles.z[i] = (float(rand())/RAND_MAX);
       //printf("Particle position: %f,%f,%f\n",particles.x[i],particles.y[i],particles.z[i]);

       particles.vx[i] = (float(rand())/RAND_MAX);
       particles.vy[i] = (float(rand())/RAND_MAX);
       particles.vz[i] = (float(rand())/RAND_MAX);
    }
}

void Basic_SPH_System::setSimState(SimState state){
    this->simState = state;
    
    this->h9 = pow(state.h, 9);
    this->h6 = pow(state.h, 6);
    this->h3 = pow(state.h, 3);
    this->h2 = pow(state.h, 2);
    
    this->h9_315 = 315./(64*M_PI*h9);
    this->h6_15  = 15./(    M_PI * h6);
    this->h3_15  = 15./(2 * M_PI * h3);
}

float Basic_SPH_System::distanceIJ(int i, int j){
    float rx = pow(this->particles.x[i] - this->particles.x[j],2);
    float ry = pow(this->particles.y[i] - this->particles.y[j],2);
    float rz = pow(this->particles.z[i] - this->particles.z[j],2);
    return sqrt(rx + ry + rz);
}

float Basic_SPH_System::evalKernel_poly6(int i, int j, float h){
    float r = this->distanceIJ(i,j);
    
    if ( r <= h){
        return h9_315 * pow( h2 - r*r, 3);
    } else {
        return 0;
    }
}

float Basic_SPH_System::evalKernel_spiky(int i, int j, float h){
    float r = this->distanceIJ(i, j);
    
    if (r <= h){
        return h6_15 * pow(h - r, 3);
    } else {
        return 0;
    }
}

float Basic_SPH_System::evalKernel_visc(int i, int j, float h){
    float r = this->distanceIJ(i, j);
    if (r <= h){
        return h3_15 * ( - r*r*r/(2*h3) + r*r/(h2) + h/(2*r) - 1);
    } else {
        return 0;
    }
}


void Basic_SPH_System::update_velocities_dummy(float dt)
{
    for( size_t i = 0; i < particles.n_particles; i++ ) {
        update_particle_velocity_dummy(i,dt);
    }

}

void Basic_SPH_System::update_positions_dummy(float dt)
{
    for( size_t i = 0; i < particles.n_particles; i++ ) {
        update_particle_position_dummy(i,dt);
    }
}

void Basic_SPH_System::calculate_Densities(){
    
    for(int i = 0; i < this->particles.n_particles; i++){
        
    }
}

void Basic_SPH_System::update_particle_position_dummy(int i, float dt)
{
    particles.x[i] = particles.x[i] + particles.vx[i]*dt;
    particles.y[i] = particles.y[i] + particles.vy[i]*dt;
    particles.z[i] = particles.z[i] + particles.vz[i]*dt;
}

//Bounce around 0-1 box
void Basic_SPH_System::update_particle_velocity_dummy(int i, float dt)
{
    if( particles.x[i] <= 0 or particles.x[i] >= 1 ) {
        particles.vx[i] = -particles.vx[i];
    }
    if( particles.y[i] <= 0 or particles.y[i] >= 1 ) {
        particles.vy[i] = -particles.vy[i];
    }
    if( particles.z[i] <= 0 or particles.z[i] >= 1 ) {
        particles.vz[i] = -particles.vz[i];
    }
    
}

size_t Basic_SPH_System::get_particle_number()
{
    return particles.n_particles;
}

//Main function that updates our particle system
void Basic_SPH_System::run_step(float dt)
{
    //Update the neighbor information
    uniform_grid.build(particles.x,particles.y,particles.z,particles.n_particles);

    //This way we can choose to have a function per particle
    //Or a general one that does all of them
    update_velocities_dummy(dt);
    update_positions_dummy(dt);
    
}


Basic_SPH_System::~Basic_SPH_System(){
    free(particles.Fx);
    free(particles.Fy);
    free(particles.Fz);
    
    free(particles.x);
    free(particles.y);
    free(particles.z);
    
    free(particles.vx);
    free(particles.vy);
    free(particles.vz);

    free(particles.rho);
}
