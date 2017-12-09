#include "Basic_SPH_System.h"
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <math.h>
#include <thinks/poissonDiskSampling.hpp>
#include "Particle_Generator.h"


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

Basic_SPH_System::Basic_SPH_System(
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

    //Initialize the liquid particles by sampling
    Vec<float,3> x_min, x_max;
    x_min[0] = 0.0f;
    x_min[1] = 0.0f;
    x_min[2] = 0.0f;

    x_max[0] = 0.2;
    x_max[1] = 0.2;
    x_max[2] = 0.2;
    uint32_t max_sample_attempts = 30;
    uint32_t seed = 1981;
    std::vector<Vec<float,3>> liquid_samples = thinks::poissonDiskSampling(this->simState.h*0.7,x_min, x_max,max_sample_attempts, seed);

    particles.n_liquid_particles = liquid_samples.size();

    //Initialize the boudaries by sampling too. Have to do it for the three planes
    float sampling_distance_boundary = 0.6;
    //XY
    Vec<float,2> b_min_xy, b_max_xy;
    b_min_xy[0] = b_min_x;
    b_min_xy[1] = b_min_y;
    b_max_xy[0] = b_max_x;
    b_max_xy[1] = b_max_y;
    std::vector<Vec<float,2>> b_samples_xy = thinks::poissonDiskSampling(this->simState.h*sampling_distance_boundary, b_min_xy, b_max_xy, max_sample_attempts, seed);

    //XZ
    Vec<float,2> b_min_xz, b_max_xz;
    b_min_xz[0] = b_min_x;
    b_min_xz[1] = b_min_z;
    b_max_xz[0] = b_max_x;
    b_max_xz[1] = b_max_z;
    std::vector<Vec<float,2>> b_samples_xz = thinks::poissonDiskSampling(this->simState.h*sampling_distance_boundary, b_min_xz, b_max_xz, max_sample_attempts, seed);

    //YZ
    Vec<float,2> b_min_yz, b_max_yz;
    b_min_yz[0] = b_min_y;
    b_min_yz[1] = b_min_z;
    b_max_yz[0] = b_max_y;
    b_max_yz[1] = b_max_z;
    std::vector<Vec<float,2>> b_samples_yz = thinks::poissonDiskSampling(this->simState.h*sampling_distance_boundary, b_min_yz, b_max_yz, max_sample_attempts, seed);

    particles.n_boundary_particles = 2*(b_samples_xy.size() + b_samples_xz.size() + b_samples_yz.size() );


    //Initialize the movable particle object
    std::vector<float> x_mob;
    std::vector<float> y_mob;
    std::vector<float> z_mob;
    size_t n_mobile_particles;

    //generate_particle_cube(0.1,0.25, x_mob, y_mob, z_mob, n_mobile_particles);
    
    load_model_data(0.05, x_mob, y_mob, z_mob, n_mobile_particles);

    size_t mobile_offset = particles.n_liquid_particles + particles.n_boundary_particles;
    particles.n_mobile_particles = n_mobile_particles;









    //Reserve memory for all the particles
    particles.n_total_particles = particles.n_liquid_particles + particles.n_boundary_particles + particles.n_mobile_particles;
    
    particles.x = (float *)malloc(particles.n_total_particles*sizeof(float));
    particles.y = (float *)malloc(particles.n_total_particles*sizeof(float));
    particles.z = (float *)malloc(particles.n_total_particles*sizeof(float));
    
    particles.vx = (float *)malloc(particles.n_total_particles*sizeof(float));
    particles.vy = (float *)malloc(particles.n_total_particles*sizeof(float));
    particles.vz = (float *)malloc(particles.n_total_particles*sizeof(float));
    
    particles.Fx = (float *)malloc(particles.n_total_particles*sizeof(float));
    particles.Fy = (float *)malloc(particles.n_total_particles*sizeof(float));
    particles.Fz = (float *)malloc(particles.n_total_particles*sizeof(float));
    
    particles.rho = (float *)malloc(particles.n_total_particles*sizeof(float));

    particles.nx = (float *)malloc(particles.n_total_particles*sizeof(float));
    particles.ny = (float *)malloc(particles.n_total_particles*sizeof(float));
    particles.nz = (float *)malloc(particles.n_total_particles*sizeof(float));
    
    particles.p = (float *)malloc(particles.n_total_particles*sizeof(float));
    

    //Initialize the liquid particles
    for( size_t i = 0; i < this->particles.n_liquid_particles; i++ ) {
        particles.x[i] = liquid_samples[i][0];
        particles.y[i] = liquid_samples[i][1];
        particles.z[i] = liquid_samples[i][2];
        //printf("Particle position: %f,%f,%f\n",particles.x[i],particles.y[i],particles.z[i]);
        
        particles.vx[i] = 0.0f; //(float(rand())/RAND_MAX);
        particles.vy[i] = 0.0f; //(float(rand())/RAND_MAX);
        particles.vz[i] = 0.0f; //(float(rand())/RAND_MAX);
    }

    //Initialize the boundary particles
    size_t offset = 0;
    for( size_t i = 0; i < b_samples_xy.size(); i++ ) { //XY plane faces
        particles.x[particles.n_liquid_particles                        + i] = b_samples_xy[i][0];
        particles.y[particles.n_liquid_particles                        + i] = b_samples_xy[i][1];
        particles.z[particles.n_liquid_particles                        + i] = b_min_z;
        
        particles.vx[particles.n_liquid_particles                       + i] = 0.0f; //(float(rand())/RAND_MAX);
        particles.vy[particles.n_liquid_particles                       + i] = 0.0f; //(float(rand())/RAND_MAX);
        particles.vz[particles.n_liquid_particles                       + i] = 0.0f; //(float(rand())/RAND_MAX);

        particles.x[particles.n_liquid_particles  + b_samples_xy.size() + i] = b_samples_xy[i][0];
        particles.y[particles.n_liquid_particles  + b_samples_xy.size() + i] = b_samples_xy[i][1];
        particles.z[particles.n_liquid_particles  + b_samples_xy.size() + i] = b_max_z;
        
        particles.vx[particles.n_liquid_particles + b_samples_xy.size() + i] = 0.0f; //(float(rand())/RAND_MAX);
        particles.vy[particles.n_liquid_particles + b_samples_xy.size() + i] = 0.0f; //(float(rand())/RAND_MAX);
        particles.vz[particles.n_liquid_particles + b_samples_xy.size() + i] = 0.0f; //(float(rand())/RAND_MAX);
    }
    offset += 2*b_samples_xy.size();

    for( size_t i = 0; i < b_samples_xz.size(); i++ ) { //XY plane faces
        particles.x[particles.n_liquid_particles                        + offset + i] = b_samples_xz[i][0];
        particles.y[particles.n_liquid_particles                        + offset + i] = b_min_y;
        particles.z[particles.n_liquid_particles                        + offset + i] = b_samples_xz[i][1];
        
        particles.vx[particles.n_liquid_particles +                     + offset + i] = 0.0f; //(float(rand())/RAND_MAX);
        particles.vy[particles.n_liquid_particles +                     + offset + i] = 0.0f; //(float(rand())/RAND_MAX);
        particles.vz[particles.n_liquid_particles +                     + offset + i] = 0.0f; //(float(rand())/RAND_MAX);

        particles.x[particles.n_liquid_particles  + b_samples_xz.size() + offset + i] = b_samples_xy[i][0];
        particles.y[particles.n_liquid_particles  + b_samples_xz.size() + offset + i] = b_max_y;
        particles.z[particles.n_liquid_particles  + b_samples_xz.size() + offset + i] = b_samples_xy[i][1];;
        
        particles.vx[particles.n_liquid_particles + b_samples_xz.size() + offset + i] = 0.0f; //(float(rand())/RAND_MAX);
        particles.vy[particles.n_liquid_particles + b_samples_xz.size() + offset + i] = 0.0f; //(float(rand())/RAND_MAX);
        particles.vz[particles.n_liquid_particles + b_samples_xz.size() + offset + i] = 0.0f; //(float(rand())/RAND_MAX);
    }
    offset += 2*b_samples_xz.size();
    for( size_t i = 0; i < b_samples_yz.size(); i++ ) { //XY plane faces
        particles.x[particles.n_liquid_particles                        + offset + i] = b_min_x;
        particles.y[particles.n_liquid_particles                        + offset + i] = b_samples_yz[i][0];
        particles.z[particles.n_liquid_particles                        + offset + i] = b_samples_yz[i][1];
        
        particles.vx[particles.n_liquid_particles +                     + offset + i] = 0.0f; //(float(rand())/RAND_MAX);
        particles.vy[particles.n_liquid_particles +                     + offset + i] = 0.0f; //(float(rand())/RAND_MAX);
        particles.vz[particles.n_liquid_particles +                     + offset + i] = 0.0f; //(float(rand())/RAND_MAX);

        particles.x[particles.n_liquid_particles  + b_samples_yz.size() + offset + i] = b_max_x;
        particles.y[particles.n_liquid_particles  + b_samples_yz.size() + offset + i] = b_samples_yz[i][0];
        particles.z[particles.n_liquid_particles  + b_samples_yz.size() + offset + i] = b_samples_yz[i][1];
        
        particles.vx[particles.n_liquid_particles + b_samples_yz.size() + offset + i] = 0.0f; //(float(rand())/RAND_MAX);
        particles.vy[particles.n_liquid_particles + b_samples_yz.size() + offset + i] = 0.0f; //(float(rand())/RAND_MAX);
        particles.vz[particles.n_liquid_particles + b_samples_yz.size() + offset + i] = 0.0f; //(float(rand())/RAND_MAX);
    }
    //Initialize the mobile particles
    //Offsets for the mobile particles position (remember that at start they are centered around zero)
    float x_offset = 0.5f;
    float y_offset = 0.2f;
    float z_offset = -0.3f;
    for( int i = 0; i < n_mobile_particles; i++ ) {
        particles.x[mobile_offset + i] = x_mob[i] + x_offset;
        particles.y[mobile_offset + i] = y_mob[i] + y_offset;
        particles.z[mobile_offset + i] = z_mob[i] + z_offset;

        particles.vx[mobile_offset + i] = 0.f;
        particles.vy[mobile_offset + i] = 0.f;
        particles.vz[mobile_offset + i] = 0.f;
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

    this->h9_32pi = 32. / (M_PI * h9);
    this->h6_64 = h6 / 64.;
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

float Basic_SPH_System::evalKernel_poly6_gradient(int &i, int &j, float &h) {
    float r = this->distanceIJ(i,j);

    if ( r <= h){
        return -h9_315 * 6 * r * pow( h2 - r*r, 2);
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


float Basic_SPH_System::evalC_spline(int &i, int &j, float &h) {
    float r = this->distanceIJ(i, j);

    if ( r <= h && r > h * 0.5f){
        return h9_32pi * powf(h*r - r*r,3);
    }
    else if ( r > 0.f && r <= h * 0.5f ){
        return h9_32pi * (2 * powf(h*r - r*r,3) - h6_64 );
    }

    return 0.f;
}

void Basic_SPH_System::update_velocities_dummy(float dt)
{
    for( size_t i = 0; i < particles.n_liquid_particles; i++ ) {
        update_particle_velocity_dummy(i,dt);
    }

}

void Basic_SPH_System::update_positions_dummy(float dt)
{
    for( size_t i = 0; i < particles.n_liquid_particles; i++ ) {
        update_particle_position_dummy(i,dt);
    }
}

void Basic_SPH_System::calculate_Densities(){

    std::vector<size_t> near_cells;
    // iterate ALL particles in the system
    for(int i = 0; i < this->particles.n_liquid_particles; i++){
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

void Basic_SPH_System::calculate_Curvatures() {
    std::vector<size_t> near_cells;
    // iterate ALL particles in the system
    for(int i = 0; i < this->particles.n_liquid_particles; i++){
        this->particles.nx[i] = 0;
        this->particles.ny[i] = 0;
        this->particles.nz[i] = 0;

        near_cells.clear();
        this->uniform_grid.query_neighbors(
                this->particles.x[i],
                this->particles.y[i],
                this->particles.z[i],
                near_cells
        );

        // iterate cells

        for (auto cellId : near_cells){
            size_t *cell     = this->uniform_grid.cells[    cellId];
            size_t cell_size = this->uniform_grid.cell_size[cellId];

            Vector3T<float> ijVec;
            for (int j = 0; j < cell_size; j++){
                auto particleID = static_cast<int>(cell[j]); // in equations that's our J index
                ijVec = this->ij_vector(i, particleID);

                float n_val = this->particles.mass / this->particles.rho[particleID] * this->evalKernel_poly6_gradient(i,particleID,simState.h);
                this->particles.nx[i] += n_val * ijVec.x();
                this->particles.ny[i] += n_val * ijVec.y();
                this->particles.nz[i] += n_val * ijVec.z();
            }
        }

        this->particles.nx[i] *= this->simState.h;
        this->particles.ny[i] *= this->simState.h;
        this->particles.nz[i] *= this->simState.h;
    }

}

void Basic_SPH_System::calculate_Pressures() {
    float rho0 = this->simState.rho0;
    float k = this->simState.k;
    
    for(int i = 0; i < this->particles.n_liquid_particles; i++){
        float rho = this->particles.rho[i];
        this->particles.p[i] = (float) std::max( (float)(k * (pow(rho / rho0, 7) - 1) ), 0.0f);
//        this->particles.p[i] = (float) std::max( (float)(k * (rho - rho0)), 0.0f);
    }
}

void Basic_SPH_System::calculate_Forces() {

    std::vector<size_t> near_cells;
    // iterate ALL particles in the system
    for(int i = 0; i < this->particles.n_liquid_particles; i++){

        this->particles.Fx[i] = 0.f;
        this->particles.Fy[i] = 0.f;
        this->particles.Fz[i] = 0.f;

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

        float F_coh_x = 0.0f;
        float F_coh_y = 0.0f;
        float F_coh_z = 0.0f;

        float F_cur_x = 0.0f;
        float F_cur_y = 0.0f;
        float F_cur_z = 0.0f;

        float F_st_x = 0.0f;
        float F_st_y = 0.0f;
        float F_st_z = 0.0f;

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


                ijVec.normalize();
                float F_coh = -this->gamma * mi * mj * this->evalC_spline(i,particleID, this->simState.h);
                F_coh_x = F_coh * ijVec.x();
                F_coh_y = F_coh * ijVec.y();
                F_coh_z = F_coh * ijVec.z();

                F_cur_x = -this->gamma * mi * ( this->particles.nx[i] - this->particles.nx[particleID] );
                F_cur_y = -this->gamma * mi * ( this->particles.ny[i] - this->particles.ny[particleID] );
                F_cur_z = -this->gamma * mi * ( this->particles.nz[i] - this->particles.nz[particleID] );


                float Kij = 2 * this->simState.rho0 / (particles.rho[i] + particles.rho[particleID]);
                F_st_x += Kij * (F_coh_x + F_cur_x);
                F_st_y += Kij * (F_coh_y + F_cur_y);
                F_st_z += Kij * (F_coh_z + F_cur_z);

                // Calc Force from viscosity
                float F_visc = this->simState.mu *this->evalKernel_visc_laplacian(i,particleID,this->simState.h);

                F_visc_x +=  mi * F_visc / rhoi * (this->particles.vx[particleID] - this->particles.vx[i]) / rhoj;
                F_visc_y +=  mi * F_visc / rhoi * (this->particles.vy[particleID] - this->particles.vy[i]) / rhoj;
                F_visc_z +=  mi * F_visc / rhoi * (this->particles.vz[particleID] - this->particles.vz[i]) / rhoj;
            }
        }
        
        this->particles.Fx[i] = F_press_x  +  F_visc_x + F_st_x;
        this->particles.Fy[i] = F_press_y  +  F_visc_y + F_st_y + F_ext_y;
        this->particles.Fz[i] = F_press_z  +  F_visc_z + F_st_z;

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
    return particles.n_liquid_particles;
}

//Main function that updates our particle system
void Basic_SPH_System::run_step(float dt)
{

    move_solid_object( 0., 0., 1*dt);

    //Update the neighbor information
    uniform_grid.build(particles.x,particles.y,particles.z,particles.n_liquid_particles);
    //uniform_grid.build(particles.x,particles.y,particles.z,particles.n_total_particles);

    //This way we can choose to have a function per particle
    //Or a general one that does all of them
//    update_velocities_dummy(dt);
//    update_positions_dummy(dt);

    this->calculate_Densities();
    this->calculate_Pressures();
    this->calculate_Curvatures();
    this->calculate_Forces();
    
    float m  = this->particles.mass;
    
    for(int i = 0; i < this->particles.n_liquid_particles; i++){
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
        if (this->particles.x[i] < 0 ) {        particles.x[i] = 0.0001f; particles.vx[i] = 0.0f; }
        if (this->particles.x[i] > b_max_x )  { particles.x[i] = b_max_x; particles.vx[i]  = 0.0f; }

        if (this->particles.y[i] < 0 ) {       particles.y[i] = 0.0001f;     particles.vy[i] = 0.; }
        if (this->particles.y[i] > b_max_y)  { particles.y[i] = b_max_y;  particles.vy[i] = 0.; }

        if (this->particles.z[i] < 0 ) {       particles.z[i] = 0.0001f; particles.vz[i] = 0.0f; }
        if (this->particles.z[i] > b_max_z)  { particles.z[i] = b_max_z;  particles.vz[i] = 0.0f; }
    }
}

void Basic_SPH_System::move_solid_object( float x, float y, float z )
{
    for( size_t i = particles.n_liquid_particles + particles.n_boundary_particles;
            i < particles.n_total_particles; i++ ) {
        particles.x[i] += x;
        particles.y[i] += y;
        particles.z[i] += z;
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

    free(particles.nx);
    free(particles.ny);
    free(particles.nz);

    free(particles.rho);
    free(particles.p);
}
