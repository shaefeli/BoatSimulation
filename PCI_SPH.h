//
// Created by Sergey Ivannikov on 10/12/2017.
//

#ifndef PBS_PROJECT_PCI_SPH_H
#define PBS_PROJECT_PCI_SPH_H


#include "Uniform_Grid.h"
#include "Vector3T.h"
#include "SupportingStructures.h"
#include "OpenGL_Renderer.h"

class PCI_SPH {

private:

    vector<Vec<float, 3>> p_createSamplingLiquid() const;
    vector<Vec<float, 3>> p_createCubicLiquid() const;

    void p_initializeBoundary();
    void implyBCOnParticle(int i);

//    void calculate_Forces();
//    void calculate_Pressures();
//    void calculate_Curvatures(); // for F _ surface tension calculation
    void calculate_Densities();
    void calculate_Forces_StageOne(); // calculate Viscouse, G, Ext, SurfaceTension forces as Stage 1
    float calculate_DensityPressureCorrection(); // returns the rho_err
    void calculate_Forces_Pressure();
    /**
     Different kernels for different properties
     */
    float evalKernel_poly6(int &i, int &j, float &h);
    float evalKernel_poly6_gradient(int &i, int &j, float &h);

    float evalKernel_spiky(int i, int j, float h);
    float evalkernel_spiky_gradient(int &i, int &j, float &h);

    float evalKernel_visc( int i, int j, float h);
    float evalKernel_visc_laplacian(int &i, int &j, float &h);

    float evalC_spline(int &i, int &j, float &h);

    // distance between particles
    Vector3T<float> ij_vectorLocation(int i, int j);
    Vector3T<float> ij_vectorVelocity(int i, int j);
    float distanceIJ(int i, int j);

    float h9,h6, h3, h2;
    float factor_poly6; // constants for kernels
    float factor_spiky;
    float factor_spiky_grad;
    float h3_15;
    float factor_visc_lapl;
    float h9_32pi;
    float h6_64;

    float gamma = 0.5; // Surface tension coefficient
    /**
     * // Delta from PCI SPH paper
     *              -1
     * = ----------------------------------
     *   betta * (-summ(grad)*summ(grad) - summ(grad * grad) )
     *
     *  where grad = gradient of the spiky kernel
     */
    float DELTA = 0.0f;

    float  current_time = 0.0f;
    size_t current_iteration = 0;


    BoundaryBox                     bBox;
    ParticlesInitialSpawningBox     iBox;
    SimState                        simState;
    UniformGridSplit                uGridSplit;


public:
    //TODO: remove this sht
    OpenGL_Renderer* debugRender;
    void debugDeltaValue();
    void precalculateDeltaValue();

    ParticlesSystemData             particles;
    std::shared_ptr<Uniform_Grid>   uniform_grid;

    PCI_SPH( BoundaryBox bBox,
             ParticlesInitialSpawningBox iBox,
             SimState simState,
             UniformGridSplit uGridSplit
    );
    ~PCI_SPH();

    //run one step of the simulation ( could adjust time step )
    void run_step();


    void preCalculateCoefficients();
    void updateParticlesPositionAndVelocity();

    float getCurrentTime() const { return this->current_time; };
    size_t getCurrentIteration() const { return this->current_iteration;};
};


#endif //PBS_PROJECT_PCI_SPH_H
