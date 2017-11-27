#include <cstddef>

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
    public:
        particle_information_t particles;
        Basic_SPH_System( size_t n_particles );
        void run_step(float dt);

};
