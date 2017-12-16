//
// Created by Sergey Ivannikov on 10/12/2017.
//

#include "PCI_SPH.h"
#include <thinks/poissonDiskSampling.hpp>
#include "SupportingStructures.h"
#include "Particle_Generator.h"
#include <algorithm>    // std::max

PCI_SPH::PCI_SPH(
        BoundaryBox                 bBox,
        ParticlesInitialSpawningBox iBox,
        SimState                    simState,
        UniformGridSplit            uGridSplit) :
        bBox(bBox), iBox(iBox), simState(simState), uGridSplit(uGridSplit)
{
    this->particles.mass = simState.mass;
    
    /*
    this->uniform_grid = std::make_shared<Uniform_Grid>(
            bBox.x1, bBox.y1, bBox.z1,
            bBox.x2, bBox.y2, bBox.z2,
            uGridSplit.cells_x,uGridSplit.cells_y,uGridSplit.cells_z);
            */
    //std::cerr<<"initializing uniform grid"<<std::endl;
    this->uniform_grid = new Uniform_Grid(
            bBox.x1, bBox.y1, bBox.z1,
            bBox.x2, bBox.y2, bBox.z2,
            uGridSplit.cells_x,uGridSplit.cells_y,uGridSplit.cells_z);
    //std::cerr<<"end initializing uniform grid"<<std::endl;


    /**
     *  Sampling for liquid initial positions
     */

    //TODO: LOOK HERE
    vector<Vec<float,3>> liquid_samples =  p_createCubicLiquid();
//    vector<Vec<float,3>> liquid_samples = p_createSamplingLiquid();
//    vector<Vec<float,3>> liquid_samples;
//    Vec<float,3> p1,p2,p3;
//    p1[0] = 0.1f;  p1[1] = 0.0f;   p1[2] = 0.1f;
//    p2[0] = 0.25f; p2[1] = 0.001f; p2[2] = 0.1f;
//    p3[0] = 0.4f;  p3[1] = 0.001f; p3[2] = 0.1f;
//    liquid_samples.push_back(p1);
//    liquid_samples.push_back(p2);
//    liquid_samples.push_back(p3);
    this->particles.n_liquid_particles_start = 0;
    this->particles.n_liquid_particles       = static_cast<unsigned int>(liquid_samples.size());


    /**
     *  Initialize boundary particles (might be empty for now)
     */
    p_initializeBoundary(); // for now ( and maybe forever) we treat boundary box by


    /**
     *  Load BOAT info - do initialize it here (see below)
     */
    std::vector<float> x_mob;
    std::vector<float> y_mob;
    std::vector<float> z_mob;
    size_t n_mobile_particles = 0;

    double scale = 0.1;
    load_model_data(0.03, scale, x_mob, y_mob, z_mob, n_mobile_particles);
    //this->generate_particle_cube(.2f, 0.02, x_mob, y_mob,z_mob,n_mobile_particles);
    int mobile_offset = particles.n_liquid_particles + particles.n_boundary_particles;
    this->particles.n_mobile_particles_start = particles.n_boundary_particles_start + particles.n_boundary_particles;
    this->particles.n_mobile_particles = n_mobile_particles;

    /**
     *  Total number of particles now is known - we can allocate and fill all particles info
     */
    this->particles.n_total_particles = particles.n_liquid_particles   +
                                  particles.n_boundary_particles +
                                  particles.n_mobile_particles;
    //std::cout<<"Allocating for n particles:"<<(this->particles.n_total_particles)<<std::endl;

    this->particles.allocateMemoryForNParticles(particles.n_total_particles);



    /**
     *  Initialize the liquid particles
     */

    for( size_t i = this->particles.n_liquid_particles_start;
         i < this->particles.n_liquid_particles + this->particles.n_liquid_particles_start;
         i++ )
    {
        particles.x[i] = liquid_samples[i][0];
        particles.y[i] = liquid_samples[i][1];
        particles.z[i] = liquid_samples[i][2];

        particles.vx[i] = 0.0f;
        particles.vy[i] = 0.0f;
        particles.vz[i] = 0.0f;
    }


    /**
     *  Initialize the BOAT particles
     */
    float x_offset = 0.2f;
    float y_offset = 0.2f;
    float z_offset = 0.5f;
    for( size_t i = particles.n_mobile_particles_start;
                i < particles.n_mobile_particles_start + particles.n_mobile_particles;
                i++ )
    {
        
        particles.x[i] = x_mob[i-particles.n_mobile_particles_start] + x_offset;
        particles.y[i] = y_mob[i-particles.n_mobile_particles_start] + y_offset;
        particles.z[i] = z_mob[i-particles.n_mobile_particles_start] + z_offset;
        //std::cout<<"i:"<<i<<" "<<particles.x[i]<<","<<particles.y[i]<<","<<particles.z[i]<<std::endl;

        particles.x_star[i] = x_mob[i-particles.n_mobile_particles_start] + x_offset;
        particles.y_star[i] = y_mob[i-particles.n_mobile_particles_start] + y_offset;
        particles.z_star[i] = z_mob[i-particles.n_mobile_particles_start] + z_offset;

        particles.vx[i] = 0.f;
        particles.vy[i] = 0.f;
        particles.vz[i] = 0.f;

        particles.vx_star[i] = 0.f;
        particles.vy_star[i] = 0.f;
        particles.vz_star[i] = 0.f;

        particles.p[i] = 0.f;
        particles.p[i] = 0.f;
        particles.p[i] = 0.f;

        particles.rho[i] = simState.rho0;
    }

    /**
     * Pre Calculate all coefficients needed for Kernels computations
     */
    preCalculateCoefficients();


    /**
     *  We need to calculate densities to start PCI SPH algorithm
     */
//    this->calculate_Densities();

}

void PCI_SPH::preCalculateCoefficients() {
    h9 = powf(this->simState.kernel_radius, 9);
    h6 = powf(this->simState.kernel_radius, 6);
    h3 = powf(this->simState.kernel_radius, 3);
    h2 = powf(this->simState.kernel_radius, 2);


    factor_poly6      = static_cast<float>(315.f / (64.f * M_PI * h9));
    factor_spiky       = static_cast<float>(15.f  / (M_PI * h6));
    factor_spiky_grad  = static_cast<float>(-45.f / (M_PI * h6));

    h3_15       = static_cast<float>(15.f / (2.f * M_PI * h3));
    factor_visc_lapl  = static_cast<float>(45.f / (M_PI * h6));

    h9_32pi = static_cast<float>(32.f / (M_PI * h9));
    h6_64   = h6 / 64.f;
}

void PCI_SPH::p_initializeBoundary()  {
    this->particles.n_boundary_particles_start = particles.n_liquid_particles_start + particles.n_liquid_particles;
    this->particles.n_boundary_particles = 0;
}

vector<Vec<float, 3>> PCI_SPH::p_createSamplingLiquid() const {
    Vec<float,3> x_min, x_max;
    x_min[0] = this->iBox.x1;
    x_min[1] = this->iBox.y1;
    x_min[2] = this->iBox.z1;

    x_max[0] = this->iBox.x2;
    x_max[1] = this->iBox.y2;
    x_max[2] = this->iBox.z2;
    uint32_t max_sample_attempts = 30;
    uint32_t seed = 1981;
    vector<Vec<float,3>> liquid_samples = thinks::poissonDiskSampling(
            this->iBox.spawningRadius,
            x_min, x_max,
            max_sample_attempts,
            seed);
    return liquid_samples;
}

PCI_SPH::~PCI_SPH() {

}

float PCI_SPH::distanceIJ(int i, int j) {
    return hypot(hypot(this->particles.x[i] - this->particles.x[j],
                       this->particles.y[i] - this->particles.y[j]),
                       this->particles.z[i] - this->particles.z[j]);
}

Vector3T<float> PCI_SPH::ij_vectorLocation(int i, int j) {
    return Vector3T<float>(
            this->particles.x[i] - this->particles.x[j],
            this->particles.y[i] - this->particles.y[j],
            this->particles.z[i] - this->particles.z[j]
    );
}

Vector3T<float> PCI_SPH::ij_vectorVelocity(int i, int j) {
    return Vector3T<float>(
            this->particles.vx[i] - this->particles.vx[j],
            this->particles.vy[i] - this->particles.vy[j],
            this->particles.vz[i] - this->particles.vz[j]
    );
}

float PCI_SPH:: evalKernel_poly6(int &i, int &j, float &h) {
    float r = this->distanceIJ(i,j);
    if ( r <= h){
        float d = h2 - r*r;
        return factor_poly6 * d*d*d;
    } else {
        return 0;
    }
}

float PCI_SPH::evalKernel_poly6_gradient(int &i, int &j, float &h) {
    float r = this->distanceIJ(i,j);
    if ( r <= h){
        float d = h2 - r*r;
        return -factor_poly6 * 6.0f * r * d*d;
    } else {
        return 0;
    }
}

float PCI_SPH::evalKernel_spiky(int i, int j, float h) {
    float r = this->distanceIJ(i, j);
    if (r <= h){
        float d = h - r;
        return factor_spiky * d*d*d;
    } else {
        return 0;
    }
}

float PCI_SPH::evalkernel_spiky_gradient(int &i, int &j, float &h) {
    float r = this->distanceIJ(i, j);
    if (r <= h){
        float d = h - r;
        return factor_spiky_grad * d*d;
    } else {
        return 0;
    }
}

float PCI_SPH::evalKernel_visc(int i, int j, float h) {
    float r = this->distanceIJ(i, j);
    if (r <= h){
        return h3_15 * ( - r*r*r/(2*h3) + r*r/(h2) + h/(2*r) - 1);
    } else {
        return 0;
    }
}

float PCI_SPH::evalKernel_visc_laplacian(int &i, int &j, float &h) {
    float r = this->distanceIJ(i, j);
    if (r <= h){
        return factor_visc_lapl * (h - r);
    } else {
        return 0;
    }
}

float PCI_SPH::evalC_spline(int &i, int &j, float &h) {
    float r = this->distanceIJ(i, j);
    if ( r <= h && r > h * 0.5f){
        float d = h*r - r*r;
        return h9_32pi * d*d*d;
    }
    else if ( r > 0.f && r <= h * 0.5f ){
        float d = h*r - r*r;
        return h9_32pi * (2.0f * d*d*d - h6_64 );
    }
    return 0.f;
}


void PCI_SPH::move_solid_object( float x, float y, float z )
{
    for( size_t i = particles.n_liquid_particles + particles.n_boundary_particles;
            i < particles.n_total_particles; i++ ) {
        particles.x[i] += x;
        particles.y[i] += y;
        particles.z[i] += z;

        mobile_mass_center_x += x;
        mobile_mass_center_y += y;
        mobile_mass_center_z += z;
    }
}


void PCI_SPH::run_step() {

    //move_solid_object( 0., 0., 1*0.002);

    //std::cerr<<"run_step:"<<std::endl;
    uniform_grid->build(particles.x, particles.y, particles.z,
                        particles.n_total_particles);
    //std::cerr<<"built"<<std::endl;


#if 1
    std::vector<size_t> near_cells;
    // Calculate RHO using position and velocity on step T
    for(int i = particles.n_liquid_particles_start; i < this->particles.n_liquid_particles; i++){
        near_cells.clear();
        this->uniform_grid->query_neighbors(
                this->particles.x[i],
                this->particles.y[i],
                this->particles.z[i],
                near_cells
        );

        // iterate cells
        int n_neighbours_active = 0;
        float rho = this->evalKernel_poly6(i,i,this->simState.kernel_radius);
        for (auto cellId : near_cells){
            size_t *cell     = this->uniform_grid->cells[    cellId];
            size_t cell_size = this->uniform_grid->cell_size[cellId];

            for (int j = 0; j < cell_size; j++){
                auto particleID = static_cast<int>(cell[j]); // in equations that's our J index
                if ( particleID == i ) continue;
                rho += this->evalKernel_poly6(i,particleID, this->simState.kernel_radius);
                n_neighbours_active++;
            }
        }
        this->particles.rho[i] = this->particles.mass * rho;
    }

    // F_visc and F_ext calculation and advections using RHO from step T
    for(int i = particles.n_liquid_particles_start; i < this->particles.n_liquid_particles; i++){

        // Force acting on the particle
        this->particles.Fx[i] = 0.f;
        this->particles.Fy[i] = 0.f;
        this->particles.Fz[i] = 0.f;

        // find neighbours
        near_cells.clear();
        this->uniform_grid->query_neighbors(
                this->particles.x[i],
                this->particles.y[i],
                this->particles.z[i],
                near_cells
        );


        /**
         * Temporary variables for computations
         */

        float F_ext_y = - this->particles.mass * this->simState.g; // GRAVITY
        float F_visc_x = 0.0f;
        float F_visc_y = 0.0f;
        float F_visc_z = 0.0f;
        float rhoi = this->particles.rho[i];

        // iterate over near cells
        for (auto cellId : near_cells){
            size_t *cell     = this->uniform_grid->cells[    cellId];
            size_t cell_size = this->uniform_grid->cell_size[cellId];

            Vector3T<float> ijVec;
            Vector3T<float> v_ijVec;
            // iterate over ELEMENTS IN THE CELL - so it is not Particle ID yet
            for (int j = 0; j < cell_size; j++){
                // get particleID to work with
                auto particleID = static_cast<int>(cell[j]); // in equations that's our J index
                if (particleID == i ) continue;

                float rhoj = this->particles.rho[particleID];
                float mj   = this->particles.mass;
                float mi   = mj;
                ijVec      = this->ij_vectorLocation(i, particleID);

                // Calc Force from viscosity
                v_ijVec = this->ij_vectorVelocity(particleID,i);
                float F_visc = this->evalKernel_visc_laplacian(i,particleID,this->simState.kernel_radius);
                F_visc   *=  mj * this->simState.mu / ( rhoj);

                F_visc_x +=   F_visc * v_ijVec.x();
                F_visc_y +=   F_visc * v_ijVec.y();
                F_visc_z +=   F_visc * v_ijVec.z();
            }
        }


        this->particles.Fx[i] = F_visc_x;
        this->particles.Fy[i] = F_visc_y + F_ext_y;
        this->particles.Fz[i] = F_visc_z;

        this->particles.vx_star[i] = particles.vx[i] + simState.dt / simState.mass * ( particles.Fx[i] );
        this->particles.vy_star[i] = particles.vy[i] + simState.dt / simState.mass * ( particles.Fy[i] );
        this->particles.vz_star[i] = particles.vz[i] + simState.dt / simState.mass * ( particles.Fz[i] );

        this->particles.x_star[i] = particles.x[i] + simState.dt * particles.vx_star[i];
        this->particles.y_star[i] = particles.y[i] + simState.dt * particles.vy_star[i];
        this->particles.z_star[i] = particles.z[i] + simState.dt * particles.vz_star[i];


        if (        particles.x_star[i] <= bBox.x1)  { particles.x_star[i] = bBox.x1; particles.vx_star[i] *= -0.5f; }
        else if (   particles.x_star[i] >= bBox.x2)  { particles.x_star[i] = bBox.x2; particles.vx_star[i] *= -0.5f; }

        if (        particles.y_star[i] <= bBox.y1)  { particles.y_star[i] = bBox.y1; particles.vy_star[i] *= -0.5f; }
        else if (   particles.y_star[i] >= bBox.y2)  { particles.y_star[i] = bBox.y2; particles.vy_star[i] *= -0.5f; }

        if (        particles.z_star[i] <= bBox.z1)  { particles.z_star[i] = bBox.z1; particles.vz_star[i] *= -0.5f; }
        else if (   particles.z_star[i] >= bBox.z2)  { particles.z_star[i] = bBox.z2; particles.vz_star[i] *= -0.5f; }
    }

    // Start to iterate with pressure to minimise error of RHO --> should be close to RHO_0 = reference density
    float rho_err = 10.f;
    float max_rho = 0.f;
    int   corr_it = 0;
    while (rho_err > 0.01 && corr_it < 10) {
        corr_it++;
        // Compute rho_star and P
        for(int i = particles.n_liquid_particles_start; i < this->particles.n_liquid_particles; i++){
            near_cells.clear();
            this->uniform_grid->query_neighbors(
                    this->particles.x_star[i],
                    this->particles.y_star[i],
                    this->particles.z_star[i],
                    near_cells
            );

            float rho = this->evalKernel_poly6(i,i,this->simState.kernel_radius);
            for (auto cellId : near_cells){
                size_t *cell     = this->uniform_grid->cells[    cellId];
                size_t cell_size = this->uniform_grid->cell_size[cellId];

                for (int j = 0; j < cell_size; j++){
                    auto particleID = static_cast<int>(cell[j]); // in equations that's our J index
                    if ( particleID == i ) continue;

                    float r =
                            pow(this->particles.x_star[i] - this->particles.x_star[particleID],2.f) +
                            pow(this->particles.y_star[i] - this->particles.y_star[particleID],2.f) +
                            pow(this->particles.z_star[i] - this->particles.z_star[particleID],2.f)
                   ;

                    float d =h2 - r;
                    if ( d >= 0){
                        rho += factor_poly6 *d*d*d;
                    }
                }
            }
            this->particles.rho[i] = this->particles.mass * rho;

            if (particles.rho[i] > max_rho ) { max_rho = particles.rho[i];}

            this->particles.p[i] = simState.k * ( pow(particles.rho[i]/simState.rho0,7.0f) - 1) ;
        }


        // Compute Pressure Force
        for(int i = particles.n_liquid_particles_start; i < this->particles.n_liquid_particles; i++){
            near_cells.clear();
            this->uniform_grid->query_neighbors(
                    this->particles.x_star[i],
                    this->particles.y_star[i],
                    this->particles.z_star[i],
                    near_cells
            );

            float rhoi = particles.rho[i];

            this->particles.Fx_p[i] = 0.0f;
            this->particles.Fy_p[i] = 0.0f;
            this->particles.Fz_p[i] = 0.0f;

            for (auto cellId : near_cells) {
                size_t *cell = this->uniform_grid->cells[cellId];
                size_t cell_size = this->uniform_grid->cell_size[cellId];

                for (int j = 0; j < cell_size; j++) {
                    auto particleID = static_cast<int>(cell[j]); // in equations that's our J index
                    if (particleID == i) continue;

                    float r =
                            pow(this->particles.x_star[i] - this->particles.x_star[particleID],2.f) +
                            pow(this->particles.y_star[i] - this->particles.y_star[particleID],2.f) +
                            pow(this->particles.z_star[i] - this->particles.z_star[particleID],2.f)
                    ;
                    r = sqrt(r);

                    if ( r <= simState.kernel_radius){
                        float rhoj = particles.rho[particleID];
//                        float Fp_scalar = simState.mass / (2 *particles.rho[j]) * (particles.p[i] + particles.p[particleID]);
                        float Fp_scalar = simState.mass * rhoi * ( particles.p[i]/(rhoi*rhoi) + particles.p[particleID]/(rhoj*rhoj)) ;
                        float d = simState.kernel_radius - r;
                        Fp_scalar *=  factor_spiky_grad * d*d;
                        Fp_scalar *= - simState.mass / particles.rho[i];


                        this->particles.Fx_p[i] += Fp_scalar *(particles.x_star[i] - particles.x_star[particleID]);
                        this->particles.Fy_p[i] += Fp_scalar *(particles.y_star[i] - particles.y_star[particleID]);
                        this->particles.Fz_p[i] += Fp_scalar *(particles.z_star[i] - particles.z_star[particleID]);
                    }
                }
            }

            particles.vx_star[i] += simState.dt/ simState.mass * particles.Fx_p[i];
            particles.vy_star[i] += simState.dt/ simState.mass * particles.Fy_p[i];
            particles.vz_star[i] += simState.dt/ simState.mass * particles.Fz_p[i];

            particles.x_star[i] += simState.dt*simState.dt / simState.mass * particles.Fx_p[i];
            particles.y_star[i] += simState.dt*simState.dt / simState.mass * particles.Fy_p[i];
            particles.z_star[i] += simState.dt*simState.dt / simState.mass * particles.Fz_p[i];
        }


        rho_err = fabsf( max_rho - simState.rho0 )/ simState.rho0;
    }
    printf("\t\t RHO_ERR = %f   max_rho = %f\n",rho_err,max_rho);

    for(int i = particles.n_liquid_particles_start; i < this->particles.n_liquid_particles; i++){
        particles.vx[i] = particles.vx_star[i];
        particles.vy[i] = particles.vy_star[i];
        particles.vz[i] = particles.vz_star[i];

        particles.x[i] = particles.x_star[i];
        particles.y[i] = particles.y_star[i];
        particles.z[i] = particles.z_star[i];

        this->implyBCOnParticle(i);
    }





#endif

#if 0
    /**
     *  Stage 1 - Compute F viscouse/g/external/curvature for all particles
     */

    this->calculate_Forces_StageOne();
    float m  = this->particles.mass;

    float eps = 1; // error should be less then 1%
    float rho_err = 100;
    size_t it = 0;
    
    
    updateParticlesPositionAndVelocity();

    printf("-----------------------------------------\n");
    while ( it < 1){
        /**
         * Predict velocity and location
         */
//        updateParticlesPositionAndVelocity();

        /**
         * Inside we predict density, its variation and update pressure.
         * Afterwards we need to compute Fp (pressure force)
         */
        this->calculate_Densities(); // recomputer densities, cause we moved particles
        rho_err = this->calculate_DensityPressureCorrection();
        rho_err = rho_err/simState.rho0; // transform into percentage difference
//        this->calculate_Forces_Pressure(); // this pressure is used on the next iteration!

        float dt = simState.dt;
        float m = particles.mass;
        
        
        for(int i = particles.n_liquid_particles_start;
            i < particles.n_liquid_particles;
            i++)
        {
            float rhoi = this->particles.rho[i];

            particles.vx[i] += dt / rhoi * ( particles.Fx_p[i]);
            particles.vy[i] += dt / rhoi * ( particles.Fy_p[i]);
            particles.vz[i] += dt / rhoi * ( particles.Fz_p[i]);
            
            particles.x[i]  += dt*dt / m * particles.Fx_p[i];
            particles.y[i]  += dt*dt / m * particles.Fy_p[i];
            particles.z[i]  += dt*dt / m * particles.Fz_p[i];
            
            implyBCOnParticle(i);
        }
        
//        this->debugRender->draw();


        printf("[%d]time  [%d]pci_it ==> rho_err = %f \n", static_cast<int>(current_iteration),it+1, rho_err);

//        printf("\t\t dR = %f \n", sqrt(pow(particles.x[0] - particles.x[1],2) + pow(particles.y[0] - particles.y[1],2) + pow(particles.z[0] - particles.z[1],2)));

        it++;
    }

    

    /**
     * Final velocity and location assignment
     */
    updateParticlesPositionAndVelocity();

    /**
     *  Update current time and iteration number
     */
#endif
    this->current_time += simState.dt;
    this->current_iteration++;
}

void PCI_SPH::updateParticlesPositionAndVelocity() {
    float dt = simState.dt;
    float m = particles.mass;
    for(int i = particles.n_liquid_particles_start;
            i < particles.n_liquid_particles;
        i++)
        {
            particles.vx[i] += dt / m * (particles.Fx[i] + particles.Fx_p[i]);
            particles.vy[i] += dt / m * (particles.Fy[i] + particles.Fy_p[i]);
            particles.vz[i] += dt / m * (particles.Fz[i] + particles.Fz_p[i]);

            particles.x[i]  += dt * particles.vx[i];
            particles.y[i]  += dt * particles.vy[i];
            particles.z[i]  += dt * particles.vz[i];

            implyBCOnParticle(i);
        }
}

void PCI_SPH::implyBCOnParticle(int i)  {
    if (        particles.x[i] <= bBox.x1)  { particles.x[i] = bBox.x1; particles.vx[i] *= 0.f; }
    else if (   particles.x[i] >= bBox.x2)  { particles.x[i] = bBox.x2; particles.vx[i] *= 0.f; }

    if (        particles.y[i] <= bBox.y1)  { particles.y[i] = bBox.y1; particles.vy[i] *= 0.f; }
    else if (   particles.y[i] >= bBox.y2)  { particles.y[i] = bBox.y2; particles.vy[i] *= 0.f; }

    if (        particles.z[i] <= bBox.z1)  { particles.z[i] = bBox.z1; particles.vz[i] *= 0.f; }
    else if (   particles.z[i] >= bBox.z2)  { particles.z[i] = bBox.z2; particles.vz[i] *= 0.f; }
}

void PCI_SPH::calculate_Densities() {
    std::vector<size_t> near_cells;
    // iterate ALL particles in the system
    for(int i = particles.n_liquid_particles_start; i < this->particles.n_liquid_particles; i++){
        near_cells.clear();
        this->uniform_grid->query_neighbors(
                this->particles.x[i],
                this->particles.y[i],
                this->particles.z[i],
                near_cells
        );

        // iterate cells
        float rho = this->evalKernel_poly6(i,i,this->simState.kernel_radius);
        for (auto cellId : near_cells){
            size_t *cell     = this->uniform_grid->cells[    cellId];
            size_t cell_size = this->uniform_grid->cell_size[cellId];

            for (int j = 0; j < cell_size; j++){
                auto particleID = static_cast<int>(cell[j]); // in equations that's our J index
                if ( particleID == i ) continue;
                rho += this->evalKernel_poly6(i,particleID, this->simState.kernel_radius);
            }
        }
        this->particles.rho[i] = this->particles.mass * rho;
    }

}

void PCI_SPH::calculate_Forces_StageOne() {
    std::vector<size_t> near_cells;
    // iterate ALL particles in the system
    for(int i = this->particles.n_liquid_particles_start;
        i < this->particles.n_liquid_particles;
        i++)
    {

        // reset some of this particle data
        // Pressure and corrisponding force
        this->particles.p[i]    = 0.0f;
        this->particles.Fx_p[i] = 0.0f;
        this->particles.Fy_p[i] = 0.0f;
        this->particles.Fz_p[i] = 0.0f;

        // Force acting on the particle
        this->particles.Fx[i] = 0.f;
        this->particles.Fy[i] = 0.f;
        this->particles.Fz[i] = 0.f;

        // Curvature info
        this->particles.nx[i] = 0;
        this->particles.ny[i] = 0;
        this->particles.nz[i] = 0;


        // find neighbours
        near_cells.clear();
        this->uniform_grid->query_neighbors(
                this->particles.x[i],
                this->particles.y[i],
                this->particles.z[i],
                near_cells
        );

        /**
         * Temporary variables for computations
         */

        float F_ext_y = - this->particles.mass * this->simState.g; // GRAVITY
        float F_visc_x = 0.0f;
        float F_visc_y = 0.0f;
        float F_visc_z = 0.0f;
        float rhoi = this->particles.rho[i];

        // iterate over near cells
        for (auto cellId : near_cells){
            size_t *cell     = this->uniform_grid->cells[    cellId];
            size_t cell_size = this->uniform_grid->cell_size[cellId];

            Vector3T<float> ijVec;
            Vector3T<float> v_ijVec;
            // iterate over ELEMENTS IN THE CELL - so it is not Particle ID yet
            for (int j = 0; j < cell_size; j++){
                // get particleID to work with
                auto particleID = static_cast<int>(cell[j]); // in equations that's our J index
                if (particleID == i ) continue;

                float rhoj = this->particles.rho[particleID];
                float mj   = this->particles.mass;
                float mi   = mj;
                ijVec      = this->ij_vectorLocation(i, particleID);

                float n_val = this->particles.mass / this->particles.rho[particleID] * this->evalKernel_poly6_gradient(i,particleID,simState.kernel_radius);
                this->particles.nx[i] += n_val * ijVec.x();
                this->particles.ny[i] += n_val * ijVec.y();
                this->particles.nz[i] += n_val * ijVec.z();


                // Calc Force from viscosity
                v_ijVec = this->ij_vectorVelocity(particleID,i);
                float F_visc = this->evalKernel_visc_laplacian(i,particleID,this->simState.kernel_radius);
                F_visc   *=  mj * this->simState.mu / ( rhoj);
                F_visc_x +=   F_visc * v_ijVec.x();
                F_visc_y +=   F_visc * v_ijVec.y();
                F_visc_z +=   F_visc * v_ijVec.z();
            }
        }

        this->particles.nx[i] *= this->simState.kernel_radius;
        this->particles.ny[i] *= this->simState.kernel_radius;
        this->particles.nz[i] *= this->simState.kernel_radius;

        this->particles.Fx[i] = F_visc_x;
        this->particles.Fy[i] = F_visc_y + F_ext_y;
        this->particles.Fz[i] = F_visc_z;

    }

    /**
     * Surface tension block!
     *
     * We repeat iterations among all particles, because only now we have info about Curvature for each particle
     */
    for(int i = this->particles.n_liquid_particles_start;
        i < this->particles.n_liquid_particles;
        i++)
    {

        // find neighbours
        near_cells.clear();
        this->uniform_grid->query_neighbors(
                this->particles.x[i],
                this->particles.y[i],
                this->particles.z[i],
                near_cells
        );

        /**
         * Temporary variables for computations
         */
        float F_coh_x = 0.0f;
        float F_coh_y = 0.0f;
        float F_coh_z = 0.0f;

        float F_cur_x = 0.0f;
        float F_cur_y = 0.0f;
        float F_cur_z = 0.0f;

        float F_st_x = 0.0f;
        float F_st_y = 0.0f;
        float F_st_z = 0.0f;

        float rhoi = this->particles.rho[i];

        // iterate over near cells
        for (auto cellId : near_cells) {
            size_t *cell     = this->uniform_grid->cells[    cellId];
            size_t cell_size = this->uniform_grid->cell_size[cellId];

            Vector3T<float> ijVec;
            // iterate over ELEMENTS IN THE CELL - so it is not Particle ID yet
            for (int j = 0; j < cell_size; j++) {
                // get particleID to work with
                
                auto particleID = static_cast<int>(cell[j]); // in equations that's our J index
                if (particleID == i ) continue;
                
                float rhoj = this->particles.rho[particleID];
                float mj   = this->particles.mass;
                float mi   = mj;
                ijVec      = this->ij_vectorLocation(i, particleID);

                ijVec.normalize();
                float F_coh = -this->gamma * mi * mj * this->evalC_spline(i,particleID, this->simState.kernel_radius);
                F_coh_x = F_coh * ijVec.x();
                F_coh_y = F_coh * ijVec.y();
                F_coh_z = F_coh * ijVec.z();

                F_cur_x = -this->gamma * mi * ( this->particles.nx[i] - this->particles.nx[particleID] );
                F_cur_y = -this->gamma * mi * ( this->particles.ny[i] - this->particles.ny[particleID] );
                F_cur_z = -this->gamma * mi * ( this->particles.nz[i] - this->particles.nz[particleID] );


                float Kij = 2 * this->simState.rho0 / (rhoi + rhoj);
                F_st_x += Kij * (F_coh_x + F_cur_x);
                F_st_y += Kij * (F_coh_y + F_cur_y);
                F_st_z += Kij * (F_coh_z + F_cur_z);
            }
        }

//        this->particles.Fx[i] += F_st_x;
//        this->particles.Fy[i] += F_st_y;
//        this->particles.Fz[i] += F_st_z;
    }
}

float PCI_SPH::calculate_DensityPressureCorrection() {

    std::vector<size_t> near_cells;
//    float betta = pow(simState.dt,2.0f)  *  pow(particles.mass,2.0f) * 2.0f / pow(simState.rho0,2.0f);
    float rho_error  = 0;

    for(int i = this->particles.n_liquid_particles_start;
            i < this->particles.n_liquid_particles;
            i++)
    {

        float rhoi     = this->particles.rho[i];
        float rho_star = rhoi - simState.rho0;

        float p_corr = 0;
        p_corr       =  this->DELTA * rho_star;

        rho_error = std::max<float>(rho_star, rho_error);
        this->particles.p[i] += p_corr;
    }

    float maxRho = *std::max_element(particles.rho, particles.rho + particles.n_liquid_particles);
    float minRho = *std::min_element(particles.rho, particles.rho + particles.n_liquid_particles);

    printf("\t\t Min:Max density -> [%f,%f]\n",minRho,maxRho);
    /*
     *  Computer new Pressure Force acting on the particle
     */
    for(int i = this->particles.n_liquid_particles_start;
        i < this->particles.n_liquid_particles;
        i++)
    {

        // find neighbours
        near_cells.clear();
        this->uniform_grid->query_neighbors(
                this->particles.x[i],
                this->particles.y[i],
                this->particles.z[i],
                near_cells
        );
        /**
         * Temporary variables for computations
         */
        float rhoi = this->particles.rho[i];

        float F_press_x = 0.0f;
        float F_press_y = 0.0f;
        float F_press_z = 0.0f;
        float pi = this->particles.p[i];

        // iterate over near cells
        for (auto cellId : near_cells) {
            size_t *cell     = this->uniform_grid->cells[    cellId];
            size_t cell_size = this->uniform_grid->cell_size[cellId];

            Vector3T<float> ijVec;
            // iterate over ELEMENTS IN THE CELL - so it is not Particle ID yet
            for (int j = 0; j < cell_size; j++) {
                // get particleID to work with
                auto particleID = static_cast<int>(cell[j]); // in equations that's our J index
                if (particleID == i ) continue;
                
                float rhoj = this->particles.rho[particleID];
                float pj   = this->particles.p[particleID];
                float mj   = this->particles.mass;
                float mi   = mj;
                ijVec      = this->ij_vectorLocation(i, particleID);
                ijVec.normalize();

                float gradP =  this->evalkernel_spiky_gradient(i, particleID, this->simState.kernel_radius);
                gradP *=   mi * mj * ( pi /(rhoi * rhoi) + pj / (rhoj * rhoj) );
//                gradP *=  gradP * 0.5 * (pi + pj) * mj / rhoj ;

                F_press_x += - gradP * ijVec.x();
                F_press_y += - gradP * ijVec.y();
                F_press_z += - gradP * ijVec.z();
            }
        }

        this->particles.Fx_p[i] = F_press_x / rhoi;
        this->particles.Fy_p[i] = F_press_y / rhoi;
        this->particles.Fz_p[i] = F_press_z / rhoi;

//        printf("\t\tFp[%d] = [%f, %f, %f]\n",i, F_press_x,F_press_y,F_press_z);
//
        if (i < 10){
            printf("\t\t\t --> dt*Fp/m = [%f, %f, %f]\n",
                   simState.dt*F_press_x/simState.mass,
                   simState.dt*F_press_y/simState.mass,
                   simState.dt*F_press_z/simState.mass);
        }
    }
    
    return rho_error;
}

void PCI_SPH::calculate_Forces_Pressure() {

}

void PCI_SPH::debugDeltaValue() {
    Vector3T<float> p1{0.f,0.f,0.f};
    Vector3T<float> p2{simState.kernel_radius * 0.7f,0.f,0.f};

    float s1 = 0.f;
    float s2 = 0.f;

    float distance = hypot(
            hypot(p1.x() - p2.x(),
                  p1.y() - p2.y()),
                  p1.z() - p2.z());

    float d = simState.kernel_radius - distance;

    float gradValue = factor_spiky_grad * d*d * distance;
}

void PCI_SPH::precalculateDeltaValue() {
    
    float poly6_d1_normalization = -945.f / (32.f * (float)M_PI * pow(simState.kernel_radius, 9.f));
    
    double value_sum[] = {0.f, 0.f, 0.f};
    double value_dot_value_sum = 0.f;
    float particle_mass = simState.rho0 / std::pow(1.f / (2.f * simState.particle_radius), 3.f);\
    float beta = pow(simState.dt,2.0f) * pow(particle_mass,2.0f) * 2 / pow(simState.rho0,2.0f);
    
    float kernel_radius = 4.f * simState.particle_radius;
    float kernel_radius2 = pow(kernel_radius,2.0f);
    float particle_size = 2.f * simState.particle_radius;
    float particle_radius = simState.particle_radius;

    float mass = 1.;
    float rho0 = 1000;
    float rhoi = 0;

    float poly6_factor = factor_poly6;
    
    for(            float z = -kernel_radius - particle_size; z <= kernel_radius + particle_size; z += particle_size) {
        for(        float y = -kernel_radius - particle_size; y <= kernel_radius + particle_size; y += particle_size) {
            for(    float x = -kernel_radius - particle_size; x <= kernel_radius + particle_size; x += particle_size) {
                // calculate poly6 diff1
                auto r = std::sqrt(x * x + y * y + z * z);
                auto r2 = r * r;
                
                if(r2 < kernel_radius2) {
                    float factor = poly6_d1_normalization * std::pow(kernel_radius2 - r2, 2.f);
                    float value[] = {
                        -factor * x,
                        -factor * y,
                        -factor * z,
                    };

                    rhoi += poly6_factor * std::pow(kernel_radius2 - r2, 3.f);
                    
                    for(int i = 0; i < 3; i++) {
                        // add value
                        value_sum[i] += value[i];
                        
                        // dot product of value
                        value_dot_value_sum += value[i] * value[i];
                    }
                    
                }
            }
        }
    }

    mass = rho0 / rhoi;

    // dot product of value sum
    float value_sum_dot_value_sum = 0.f;
    for(int i = 0; i < 3; i++) {
        value_sum_dot_value_sum += value_sum[i] * value_sum[i];
    }
    this->DELTA = -1.f / (beta * (-value_sum_dot_value_sum - value_dot_value_sum));

}

vector<Vec<float, 3>> PCI_SPH::p_createCubicLiquid() const {
    vector<Vec<float, 3>> res;
    float p_offset = iBox.spawningRadius;
    for (           float x = iBox.x1; x <= iBox.x2; x += p_offset){
        for (       float y = iBox.y1; y <= iBox.y2; y += p_offset){
            for (   float z = iBox.z1; z <= iBox.z2; z += p_offset){
                res.emplace_back(x,y,z);
            }
        }
    }
    return res;
}



void PCI_SPH::generate_particle_cube(
        float length,  //Length of cube side
        float h,        //Distance between particles
        std::vector<float> &xv,
        std::vector<float> &yv,
        std::vector<float> &zv,
        size_t &n_particles
) const
{
    float xmin = -length;
    float ymin = -length;
    float zmin = -length;
    float xmax =  length;
    float ymax =  length;
    float zmax =  length;

    size_t x_particles = ceil((xmax-xmin)/h);
    size_t y_particles = ceil((ymax-ymin)/h);
    size_t z_particles = ceil((zmax-zmin)/h);
    xv = std::vector<float>(x_particles*y_particles);
    yv = std::vector<float>(x_particles*y_particles);
    zv = std::vector<float>(x_particles*y_particles);
    int num_p = 0;
    
    for (size_t i = 0; i < x_particles; i++)
        for (size_t j = 0; j < y_particles; j++){
            float xval = i*h+xmin;
            float yval = j*h+ymin;
            xv[i*y_particles + j] = xval;
            yv[i*y_particles + j] = yval;
            zv[i*y_particles + j] = zmin;
            //std::cout<<"gen cube:"<<xval<<","<<yval<<","<<zmin<<std::endl;


//            xv[x_particles*y_particles + i*x_particles + j] = i*h+xmin;
//            yv[x_particles*y_particles + i*x_particles + j] = j*h+ymin;
//            zv[x_particles*y_particles + i*x_particles + j] = zmax;
//
//            num_p++;
            num_p++;
        }

    int start_2 = num_p;

//    for (size_t k = 0; k < z_particles; k++)
//        for (size_t j = 1; j < y_particles - 1; j++){
//            xv[start_2 + k*z_particles + j] = xmin;
//            yv[start_2 + k*z_particles + j] = j*h+ymin;
//            zv[start_2 + k*z_particles + j] = k*h+zmin;
//
//
//            xv[start_2 + z_particles*y_particles + k*z_particles + j] = xmax;
//            yv[start_2 + z_particles*y_particles + k*z_particles + j] = j*h+ymin;
//            zv[start_2 + z_particles*y_particles + k*z_particles + j] = k*h+zmin;
//
//            num_p++;
//            num_p++;
//        }
//
//    int start_3 = num_p;
//
//    for (size_t k = 1; k < z_particles - 1; k++)
//        for (size_t i = 1; i < x_particles - 1; i++){
//            xv[start_3 + k*z_particles + i] = i*h+xmin;
//            yv[start_3 + k*z_particles + i] = ymin;
//            zv[start_3 + k*z_particles + i] = k*h+zmin;
//
//
//            xv[start_3 + z_particles*x_particles + k*z_particles + i] = i*h+xmin;
//            yv[start_3 + z_particles*x_particles + k*z_particles + i] = ymax;
//            zv[start_3 + z_particles*x_particles + k*z_particles + i] = k*h+zmin;
//
//            num_p++;
//            num_p++;
//        }
    
//    for( size_t i = 0; i < x_particles; i++ ) {
//        for( size_t j = 0; j < y_particles; j++ ) {
//            for( size_t k = 0; k < z_particles; k++ ) {
//                if (i == 0 || i == x_particles-1 || j == 0 || j == y_particles-1 || k == 0 || k == z_particles-1){
//                    xv[i*y_particles*z_particles + j*z_particles + k] = i*h+xmin;
//                    yv[i*y_particles*z_particles + j*z_particles + k] = j*h+ymin;
//                    zv[i*y_particles*z_particles + j*z_particles + k] = k*h+zmin;
//                    num_p++;
//                }
//            }
//        }
//    }
    n_particles = static_cast<size_t>(num_p);
}




