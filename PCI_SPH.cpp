//
// Created by Sergey Ivannikov on 10/12/2017.
//

#include "PCI_SPH.h"
#include <thinks/poissonDiskSampling.hpp>
#include "SupportingStructures.h"
#include <algorithm>    // std::max

PCI_SPH::PCI_SPH(
        BoundaryBox                 bBox,
        ParticlesInitialSpawningBox iBox,
        SimState                    simState,
        UniformGridSplit            uGridSplit) :
        bBox(bBox), iBox(iBox), simState(simState), uGridSplit(uGridSplit)
{
    this->particles.mass = simState.mass;

    this->uniform_grid = std::make_shared<Uniform_Grid>(
            bBox.x1, bBox.y1, bBox.z1,
            bBox.x2, bBox.y2, bBox.z2,
            uGridSplit.cells_x,uGridSplit.cells_y,uGridSplit.cells_z);


    /**
     *  Sampling for liquid initial positions
     */
    vector<Vec<float,3>> liquid_samples = p_createSamplingLiquid();
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
    unsigned int n_mobile_particles = 0;
//    load_model_data(0.05, x_mob, y_mob, z_mob, n_mobile_particles);
    int mobile_offset = particles.n_liquid_particles + particles.n_boundary_particles;
    this->particles.n_mobile_particles_start = particles.n_boundary_particles_start + particles.n_boundary_particles;
    this->particles.n_mobile_particles = n_mobile_particles;

    /**
     *  Total number of particles now is known - we can allocate and fill all particles info
     */
    this->particles.n_total_particles = particles.n_liquid_particles   +
                                  particles.n_boundary_particles +
                                  particles.n_mobile_particles;

    this->particles.allocateMemoryForNParticles(particles.n_total_particles);



    /**
     *  Initialize the liquid particles
     */

    for( size_t i = this->particles.n_liquid_particles_start;
         i < this->particles.n_liquid_particles;
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
    float x_offset = 0.5f;
    float y_offset = 0.2f;
    float z_offset = -0.3f;
    for( size_t i = particles.n_mobile_particles_start;
         i < particles.n_mobile_particles_start + n_mobile_particles;
         i++ )
    {
        particles.x[i] = x_mob[i] + x_offset;
        particles.y[i] = y_mob[i] + y_offset;
        particles.z[i] = z_mob[i] + z_offset;

        particles.vx[i] = 0.f;
        particles.vy[i] = 0.f;
        particles.vz[i] = 0.f;

        particles.rho[i] = simState.rho0 * 2;
    }

    /**
     * Pre Calculate all coefficients needed for Kernels computations
     */
    preCalculateCoefficients();


    /**
     *  We need to calculate densities to start PCI SPH algorithm
     */
    this->calculate_Densities();

}

void PCI_SPH::preCalculateCoefficients() {
    h9 = powf(this->simState.kernel_radius, 9);
    h6 = powf(this->simState.kernel_radius, 6);
    h3 = powf(this->simState.kernel_radius, 3);
    h2 = powf(this->simState.kernel_radius, 2);


    h9_315      = static_cast<float>(315.f / (64.f * M_PI * h9));
    h6_15       = static_cast<float>(15.f  / (M_PI * h6));
    h6_15_grad  = static_cast<float>(-45.f / (M_PI * h6));

    h3_15       = static_cast<float>(15.f / (2.f * M_PI * h3));
    h3_15_visc  = static_cast<float>(45.f / (M_PI * h6));

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
        return h9_315 * d*d*d;
    } else {
        return 0;
    }
}

float PCI_SPH::evalKernel_poly6_gradient(int &i, int &j, float &h) {
    float r = this->distanceIJ(i,j);
    if ( r <= h){
        float d = h2 - r*r;
        return -h9_315 * 6.0f * r * d*d;
    } else {
        return 0;
    }
}

float PCI_SPH::evalKernel_spiky(int i, int j, float h) {
    float r = this->distanceIJ(i, j);
    if (r <= h){
        float d = h - r;
        return h6_15 * d*d*d;
    } else {
        return 0;
    }
}

float PCI_SPH::evalkernel_spiky_gradient(int &i, int &j, float &h) {
    float r = this->distanceIJ(i, j);
    if (r <= h){
        float d = h - r;
        return h6_15_grad * d*d;
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
        return h3_15_visc * (h - r);
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



void PCI_SPH::run_step() {


    /**
     * To boost computations we do it only once even if we have
     * use predictor-corrector-like scheme
     */
    uniform_grid->build(particles.x, particles.y, particles.z,
                        particles.n_total_particles);



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
        this->calculate_Forces_Pressure(); // this pressure is used on the next iteration!

        float dt = simState.dt;
        float m = particles.mass;
        
        
        for(int i = particles.n_liquid_particles_start;
            i < particles.n_liquid_particles;
            i++)
        {
            particles.vx[i] += dt / m * ( particles.Fx_p[i]);
            particles.vy[i] += dt / m * ( particles.Fy_p[i]);
            particles.vz[i] += dt / m * ( particles.Fz_p[i]);
            
            particles.x[i]  += dt * particles.vx[i];
            particles.y[i]  += dt * particles.vy[i];
            particles.z[i]  += dt * particles.vz[i];
            
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
//    updateParticlesPositionAndVelocity();

    /**
     *  Update current time and iteration number
     */
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
    if (        particles.x[i] < bBox.x1)  { particles.x[i] = bBox.x1 + 0.001f; particles.vx[i] *= 0.f; }
    else if (   particles.x[i] > bBox.x2)  { particles.x[i] = bBox.x2 - 0.001f; particles.vx[i] *= 0.f; }

    if (        particles.y[i] < bBox.y1)  { particles.y[i] = bBox.y1 + 0.001f; particles.vy[i] *= 0.f; }
    else if (   particles.y[i] > bBox.y2)  { particles.y[i] = bBox.y2 - 0.001f; particles.vy[i] *= 0.f; }

    if (        particles.z[i] < bBox.z1)  { particles.z[i] = bBox.z1 + 0.001f; particles.vz[i] *= 0.f; }
    else if (   particles.z[i] > bBox.z2)  { particles.z[i] = bBox.z2 - 0.001f; particles.vz[i] *= 0.f; }
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


        if (i  < 10){
            printf("\t\trho[%d] = %f\n",i, particles.rho[i]);
        }
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

//        // find neighbours
//        near_cells.clear();
//        this->uniform_grid->query_neighbors(
//                this->particles.x[i],
//                this->particles.y[i],
//                this->particles.z[i],
//                near_cells
//        );
//        /**
//         * Temporary variables for computations
//         */
//        float rhoi = this->particles.rho[i];
//
//        /**
//         *  Temp summs to hold Summ( grad_Wij ) and Summ( (grad_Wij)^2 )
//         */
//
//        float s1x = 0.0f;
//        float s1y = 0.0f;
//        float s1z = 0.0f;
//
//        float s2  = 0.0f;
//
//        // iterate over near cells
//        for (auto cellId : near_cells) {
//            size_t *cell     = this->uniform_grid->cells[    cellId];
//            size_t cell_size = this->uniform_grid->cell_size[cellId];
//
//            Vector3T<float> ijVec;
//            // iterate over ELEMENTS IN THE CELL - so it is not Particle ID yet
//            for (int j = 0; j < cell_size; j++) {
//                // get particleID to work with
//                auto particleID = static_cast<int>(cell[j]); // in equations that's our J index
//                if (particleID == i ) continue;
//
//                ijVec      = this->ij_vectorLocation(i, particleID);
//                ijVec.normalize();
//
//                //TODO: actually I think we can make it simplier!
//                float gradVal = this->evalkernel_spiky_gradient(i, particleID, simState.kernel_radius);
//                s1x += gradVal * ijVec.x();
//                s1y += gradVal * ijVec.y();
//                s1z += gradVal * ijVec.z();
//
//                s2 += gradVal * gradVal;
//            }
//        }

        float rhoi     = this->particles.rho[i];
        float rho_star = this->particles.rho[i] - simState.rho0;

        float p_corr = 0;
        p_corr   =  this->DELTA * rho_star;

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
//                gradP *=  gradP * 0.5 * (pi + pj) * mi * mj / ( rhoi * rhoj) ;

                F_press_x += - gradP * ijVec.x();
                F_press_y += - gradP * ijVec.y();
                F_press_z += - gradP * ijVec.z();
            }
        }

        this->particles.Fx_p[i] = F_press_x;
        this->particles.Fy_p[i] = F_press_y;
        this->particles.Fz_p[i] = F_press_z;

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

    float gradValue = h6_15_grad * d*d * distance;
}

void PCI_SPH::precalculateDeltaValue() {
    
    float poly6_d1_normalization = -945.f / (32.f * (float)M_PI * pow(simState.kernel_radius, 9.f));
    
    float value_sum[] = {0.f, 0.f, 0.f};
    float value_dot_value_sum = 0.f;
    float particle_mass = simState.rho0 / std::pow(1.f / (2.f * simState.particle_radius), 3.f);\
    float beta = pow(simState.dt,2.0f) * pow(particle_mass,2.0f) * 2 / pow(simState.rho0,2.0f);
    
    float kernel_radius = 4.f * simState.particle_radius;
    float kernel_radius2 = pow(kernel_radius,2.0f);
    float particle_size = 2.f * simState.particle_radius;
    float particle_radius = simState.particle_radius;

    float mass = 1.;
    float rho0 = 1000;
    float rhoi = 0;

    float poly6_factor = h9_315;
    
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
    
    
#if 0
    // we need to pre-create small volume of particles ( can use Poisson Disk sampling )
    // they're centered around one particle
    Vector3T<float> p1{0.f,0.f,0.f};
    
    
    // by the way - we can do it in 2D dimenstion as well... should not affect computations much...
    
    Vec<float,3> x_min, x_max;
    x_min[0] = 0;
    x_min[1] = 0;
    x_min[2] = 0;
    
    x_max[0] = simState.kernel_radius * 7;
    x_max[1] = simState.kernel_radius * 7;
    x_max[2] = simState.kernel_radius * 7;
    uint32_t max_sample_attempts = 100;
    uint32_t seed = 1981;
    vector<Vec<float,3>> liquid_samples = thinks::poissonDiskSampling(
                                                                      this->iBox.spawningRadius,
                                                                      x_min, x_max,
                                                                      max_sample_attempts*0.9,
                                                                      seed);
    
    Vec<float,3> theLOC;
    
    float kernel_radius = simState.kernel_radius;
    int ind = -1;
    for ( int i = 0; i < liquid_samples.size(); i++){
        auto loc = liquid_samples[i];
        if ( loc[0] >= 2.f*kernel_radius  && loc[0] <= 4.f*kernel_radius  &&
             loc[1] >= 2.f*kernel_radius  && loc[1] <= 4.f*kernel_radius  &&
             loc[2] >= 2.f*kernel_radius  && loc[2] <= 4.f*kernel_radius )
        {
            theLOC = loc;
            ind    = i;
        }
    }


    if (ind >= 0){
        Vector3T<float> ijVec;
        float s1 = 0;
        float s2 = 0;

        for (int i = 0; i < liquid_samples.size(); i++){
            if (i == ind) continue;

            ijVec = Vector3T<float>( theLOC[0] - liquid_samples[i][0],
                                     theLOC[1] - liquid_samples[i][1],
                                     theLOC[2] - liquid_samples[i][2]) ;

            float R = ijVec.length();
            ijVec.normalize();

            float kernelGrad = 0.0f;
            if (R <= simState.kernel_radius){
                float d = kernel_radius - R;
                kernelGrad = h6_15_grad * d*d;
            }

            s1 += kernelGrad;
            s2 += kernelGrad * kernelGrad;
        }

        float betta = pow(simState.dt,2.0f) * pow(simState.mass,2.0f) * 2 / pow(simState.rho0,2.0f);
        DELTA = -1.0f / ( betta *( -s1 -s2) );
    }
#endif
}

vector<Vec<float, 3>> PCI_SPH::p_createCubicLiquid() const {
    vector<Vec<float, 3>> res;
    float p_offset = this->simState.kernel_radius * 0.7f;
    for (           float x = iBox.x1; x <= iBox.x2; x += p_offset){
        for (       float y = iBox.y1; y <= iBox.y2; y += p_offset){
            for (   float z = iBox.z1; z <= iBox.z2; z += p_offset){

            }
        }
    }

    return res;
}




