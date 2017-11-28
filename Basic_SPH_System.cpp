#include "Basic_SPH_System.h"
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <math.h>
#include <thinks/poissonDiskSampling.hpp>


template <typename T, std::size_t N>
class Vec
{
public:
    typedef T value_type;
    static const std::size_t size = N;
    Vec() {}
    T& operator[](std::size_t i) { return _data[i]; }
    const T& operator[](std::size_t i) const { return _data[i]; }
private:
    T _data[N];
};

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
    
}

void Basic_SPH_System::finilizeInit() {

    Vec<float,3> x_min, x_max;
    x_min[0] = 0.0f;
    x_min[1] = 0.0f;
    x_min[2] = 0.0f;
    x_max[0] = 0.5;
    x_max[1] = 0.5;
    x_max[2] = 0.5;
    uint32_t max_sample_attempts = 50;
    uint32_t seed = 1981;
    std::vector<Vec<float,3>> samples = thinks::poissonDiskSampling(this->simState.h,x_min, x_max,max_sample_attempts, seed);
    
    particles.n_particles = samples.size();
    
    particles.x = (float *)malloc(particles.n_particles*sizeof(float));
    particles.y = (float *)malloc(particles.n_particles*sizeof(float));
    particles.z = (float *)malloc(particles.n_particles*sizeof(float));
    
    particles.vx = (float *)malloc(particles.n_particles*sizeof(float));
    particles.vy = (float *)malloc(particles.n_particles*sizeof(float));
    particles.vz = (float *)malloc(particles.n_particles*sizeof(float));
    
    particles.Fx = (float *)malloc(particles.n_particles*sizeof(float));
    particles.Fy = (float *)malloc(particles.n_particles*sizeof(float));
    particles.Fz = (float *)malloc(particles.n_particles*sizeof(float));
    
    particles.rho = (float *)malloc(particles.n_particles*sizeof(float));
    
    particles.p = (float *)malloc(particles.n_particles*sizeof(float));
    
    for( size_t i = 0; i < this->particles.n_particles; i++ ) {
        particles.x[i] = samples[i][0];
        particles.y[i] = samples[i][1];
        particles.z[i] = samples[i][2];
        //printf("Particle position: %f,%f,%f\n",particles.x[i],particles.y[i],particles.z[i]);
        
        particles.vx[i] = 0.0f; //(float(rand())/RAND_MAX);
        particles.vy[i] = 0.0f; //(float(rand())/RAND_MAX);
        particles.vz[i] = 0.0f; //(float(rand())/RAND_MAX);
    }


    this->flag_finilized = true;
}

void Basic_SPH_System::setSimState(SimState state){
    this->simState = state;
    this->particles.mass = 0.2; //static_cast<float>(4. / 3. * M_PI * pow(state.h, 3));


    /**
     * Constants pre-calculated
     */
    this->h9 = pow(state.h, 9);
    this->h6 = pow(state.h, 6);
    this->h3 = pow(state.h, 3);
    this->h2 = pow(state.h, 2);
    
    this->h9_315     = 315./(64*M_PI*h9);
    this->h6_15      = 15./(    M_PI * h6);
    this->h6_15_grad = -45. / ( M_PI * h6);

    this->h3_15      = 15./(2 * M_PI * h3);
    this->h3_15_visc = 45./(M_PI * h6);

}

float Basic_SPH_System::distanceIJ(int i, int j){
    float rx = pow(this->particles.x[i] - this->particles.x[j],2);
    float ry = pow(this->particles.y[i] - this->particles.y[j],2);
    float rz = pow(this->particles.z[i] - this->particles.z[j],2);
    return sqrt(rx + ry + rz);
}

Vector3T<float> Basic_SPH_System::ij_vector(int i, int j) {

    return Vector3T<float>(
            this->particles.x[i] - this->particles.x[j],
            this->particles.y[i] - this->particles.y[j],
            this->particles.z[i] - this->particles.z[j]
    );
}

float Basic_SPH_System::evalKernel_poly6(int &i, int &j, float &h){
    
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

/*
 * It also NEEDS to be multiplied by vector R!! R = r(i,j) vec
 */
float Basic_SPH_System::evalkernel_spiky_gradient(int &i, int &j, float &h){
    
    float r = this->distanceIJ(i, j);

    if (r <= h){
        return h6_15_grad * pow(h - r, 2) / r;
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

float Basic_SPH_System::evalKernel_visc_laplacian(int &i, int &j, float &h) {
    
    float r = this->distanceIJ(i, j);
    if (r <= h){
        return h3_15_visc * (h - r);
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

    std::vector<size_t> near_cells;
    // iterate ALL particles in the system
    for(int i = 0; i < this->particles.n_particles; i++){
        near_cells.clear();
        this->uniform_grid.query_neighbors(
                this->particles.x[i],
                this->particles.y[i],
                this->particles.z[i],
                near_cells
        );

        // iterate cells
        float rho = this->evalKernel_poly6(i,i,this->simState.h);
        for (auto cellId : near_cells){
            size_t *cell     = this->uniform_grid.cells[    cellId];
            size_t cell_size = this->uniform_grid.cell_size[cellId];

            for (int j = 0; j < cell_size; j++){
                auto particleID = static_cast<int>(cell[j]); // in equations that's our J index
                rho += this->evalKernel_poly6(i,particleID, this->simState.h);
            }
        }
        this->particles.rho[i] = this->particles.mass * rho;
    }
}

void Basic_SPH_System::calculate_Pressures() {
    float rho0 = this->simState.rho0;
    float k = this->simState.k;
    
    for(int i = 0; i < this->particles.n_particles; i++){
        float rho = this->particles.rho[i];
//        this->particles.p[i] = (float) std::max( (float)(k * (pow(rho / rho0, 7) - 1) ), 0.0f);
        this->particles.p[i] = (float) std::max( (float)(k * (rho - rho0)), 0.0f);
    }
}

void Basic_SPH_System::calculate_Forces() {

    std::vector<size_t> near_cells;
    // iterate ALL particles in the system
    for(int i = 0; i < this->particles.n_particles; i++){
        near_cells.clear();
        this->uniform_grid.query_neighbors(
                this->particles.x[i],
                this->particles.y[i],
                this->particles.z[i],
                near_cells
        );

        // Different forces acting o i-particle
        float F_visc_x = 0.0f;
        float F_visc_y = 0.0f;
        float F_visc_z = 0.0f;

        float F_press_x = 0.0f;
        float F_press_y = 0.0f;
        float F_press_z = 0.0f;

        float F_ext_y = - this->particles.mass * this->simState.g; // GRAVITY

        float pi = this->particles.p[i];
        float rhoi = this->particles.rho[i];
        // iterate cells
        for (auto cellId : near_cells){
            size_t *cell     = this->uniform_grid.cells[    cellId];
            size_t cell_size = this->uniform_grid.cell_size[cellId];

            for (int j = 0; j < cell_size; j++){
                auto particleID = static_cast<int>(cell[j]); // in equations that's our J index
                if (particleID == i) continue; // we should not compute agains ourselves
                
                float pj   = this->particles.p[particleID];
                float rhoj = this->particles.rho[particleID];
                float mj   = this->particles.mass;
                float mi   = mj;

                // Calc Force from pressure
                Vector3T<float> ijVec = this->ij_vector(i, particleID);
//                std::cout << "[ " << i << " - " << particleID << " ] = " << sqrt(ijVec.squaredLength()) << std::endl;
                float gradP =  this->evalkernel_spiky_gradient(i, particleID, this->simState.h);
                gradP *=  rhoi * mj * ( pi /(rhoi * rhoi) + pj / (rhoj * rhoj) );
//                F_p *=  mj * ( pi + pj)/(2*rhoj);
                F_press_x += - mi/rhoi * gradP * ijVec.x();
                F_press_y += - mi/rhoi * gradP * ijVec.y();
                F_press_z += - mi/rhoi * gradP * ijVec.z();


                // Calc Force from viscosity
                float F_visc = this->simState.mu *this->evalKernel_visc_laplacian(i,j,this->simState.h);

                F_visc_x +=  mi * F_visc / rhoi * (this->particles.vx[j] - this->particles.vx[i]) / rhoj;
                F_visc_y +=  mi * F_visc / rhoi * (this->particles.vy[j] - this->particles.vy[i]) / rhoj;
                F_visc_z +=  mi * F_visc / rhoi * (this->particles.vz[j] - this->particles.vz[i]) / rhoj;
            }
        }
        
        this->particles.Fx[i] = F_press_x  +  F_visc_x;
        this->particles.Fy[i] = F_press_y  +  F_visc_y + F_ext_y;
        this->particles.Fz[i] = F_press_z  +  F_visc_z;

//        / rhoi
//        / rhoi
//        / rhoi
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
//    update_velocities_dummy(dt);
//    update_positions_dummy(dt);

    this->calculate_Densities();
    this->calculate_Pressures();
    this->calculate_Forces();
    
    float m  = this->particles.mass;
    
    for(int i = 0; i < this->particles.n_particles; i++){
        this->particles.vx[i] += dt / m * this->particles.Fx[i];
        this->particles.vy[i] += dt / m * this->particles.Fy[i];
        this->particles.vz[i] += dt / m * this->particles.Fz[i];
        
        this->particles.x[i]  += dt * this->particles.vx[i];
        this->particles.y[i]  += dt * this->particles.vy[i];
        this->particles.z[i]  += dt * this->particles.vz[i];


        /** DEBUG
         * If trying to escape - return and set velocity to 0.
         * Sort of - No-Slip boundary condition .
         */
        if (this->particles.x[i] < 0 ) {        particles.x[i] = 0.0f; particles.vx[i] = 0.0f; }
        if (this->particles.x[i] > b_max_x )  { particles.x[i] = b_max_x; particles.vx[i]  = 0.0f; }

        if (this->particles.y[i] < 0 ) {       particles.y[i] = 0.0f;     particles.vy[i] = 0.; }
        if (this->particles.y[i] > b_max_y)  { particles.y[i] = b_max_y;  particles.vy[i] = 0.; }

        if (this->particles.z[i] < 0 ) {       particles.z[i] = 0.0f; particles.vz[i] = 0.0f; }
        if (this->particles.z[i] > b_max_z)  { particles.z[i] = b_max_z;  particles.vz[i] = 0.0f; }
    }
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
    free(particles.p);
}
