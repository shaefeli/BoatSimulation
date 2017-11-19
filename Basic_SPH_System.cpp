#include "Basic_SPH_System.h"
#include <cstdio>
#include <ctime>
#include <cstdlib>

Basic_SPH_System::Basic_SPH_System( size_t n_particles )
{
    particles.n_particles = n_particles;

    particles.x = (float *)malloc(n_particles*sizeof(float));
    particles.y = (float *)malloc(n_particles*sizeof(float));
    particles.z = (float *)malloc(n_particles*sizeof(float));

    particles.vx = (float *)malloc(n_particles*sizeof(float));
    particles.vy = (float *)malloc(n_particles*sizeof(float));
    particles.vz = (float *)malloc(n_particles*sizeof(float));

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
    //This way we can choose to have a function per particle
    //Or a general one that does all of them
    update_velocities_dummy(dt);
    update_positions_dummy(dt);
    
}
