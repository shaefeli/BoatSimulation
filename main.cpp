#include "OpenGL_Renderer.h"
#include <omp.h>
#include "Basic_SPH_System.h"
#include "PCI_SPH.h"
#include "SupportingStructures.h"
#include "Maya_Interface.h"
#include "unistd.h"
#include "Particle_Generator.h"
#include <iostream>
#include <chrono>
#include <thread>

#define METHOD_BASIC_SPH 0
#define METHOD_PCI_SPH   1

void initRendererWithSimInfo(const PCI_SPH &pciSph, OpenGL_Renderer &renderer);

int main(int argc, char** argv){

    
    SimState simState;
    int startIt = 0;
    int endIt = 1000;
//    simState.dt = 0.001;
    simState.g  = 9.8;
    simState.particle_radius = 0.02f ;
    simState.kernel_radius   = 4.0f * simState.particle_radius;
    simState.k  = 1e+5;
    simState.mu = 1e-6f;
    simState.rho0 = 1000.f;
    simState.gamma = 0.f;
    simState.mass = simState.rho0 / pow( 1.f / (simState.particle_radius*2.f),3.f);
//    simState.mass = 1;
    simState.dt = 0.001;
//    simState.mass = 10.f;

#if METHOD_PCI_SPH

    BoundaryBox bBox{
            .x1 = 0.f,  .x2 = 0.4f,
            .y1 = 0.f,  .y2 = 1.0f,
            .z1 = 0.f,  .z2 = 1.0f
    };

    ParticlesInitialSpawningBox iBox{
            .x1 = 0.1f,     .x2 = 0.3f,
            .y1 = 0.1f,    .y2 = 0.9f,
            .z1 = 0.1f,     .z2 = 0.5f,
            .spawningRadius = simState.kernel_radius * 0.5f
    };
    UniformGridSplit gridSplit{
            .cells_x = simState.kernel_radius*1.0f,
            .cells_y = simState.kernel_radius*1.0f,
            .cells_z = simState.kernel_radius*1.0f
    };

    omp_set_num_threads(4);


    PCI_SPH pciSph(bBox,iBox,simState,gridSplit);



    pciSph.precalculateDeltaValue();


    OpenGL_Renderer renderer;
    initRendererWithSimInfo(pciSph, renderer);
    renderer.init(argc,argv);
    pciSph.debugRender = &renderer;

    Maya_Interface maya(pciSph.particles.n_liquid_particles,endIt, startIt);

    unsigned int it=0;

    bool render = true;
    
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    
    double elapsed_time = 0;
    size_t n_fps = 0;
    double avg_fps = 0;

    while(!glfwWindowShouldClose(renderer.getWindow())){
        std::cout << "[Simulation time:" << pciSph.getCurrentTime() << "] sec\n";
        std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
        std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        std::cout << "[Frame time:" << time_span.count() << "] sec\n";
        std::cout << "[FPS:" << avg_fps << "] sec\n";
        std::cout << pciSph.mobile_mass_center_z << "\n";
        elapsed_time += time_span.count();
        n_fps++;
        if( n_fps == 10 ) {
            avg_fps = n_fps/elapsed_time;
            //std::cout << "[FPS:" << n_fps/elapsed_time << "] sec\n";
            n_fps = 0;
            elapsed_time = 0;
        }

        t1 = std::chrono::steady_clock::now();

        if(render and it >= startIt && it <= endIt){
            maya.writeToMaya(it,
                             pciSph.particles.x,
                             pciSph.particles.y,
                             pciSph.particles.z,
                             pciSph.particles.n_liquid_particles);
            maya.writeBoatInfos(pciSph.mobile_mass_center_x,
                                pciSph.mobile_mass_center_y,
                                pciSph.mobile_mass_center_z,
                                pciSph.mobile_angle_phi,
                                pciSph.mobile_angle_theta,
                                pciSph.mobile_angle_psi,
                                it);
        }

        pciSph.run_step();
        renderer.draw();
        it++;
        ///std::this_thread::sleep_for(std::chrono::milliseconds(500));
    }


#endif


#if METHOD_BASIC_SPH
    Basic_SPH_System bsph(
                          0., 0., 0.,
                          1., 1., 1.,
                          simState.kernel_radius , simState.kernel_radius , simState.kernel_radius );
    bsph.setSimState(simState);
    bsph.finilizeInit();

    OpenGL_Renderer renderer;
    Maya_Interface maya(bsph.particles.n_liquid_particles,endIt, startIt);

    //renderer.render_info.n_liquid_particles = bsph.particles.n_liquid_particles;
    renderer.render_info.n_liquid_particles = bsph.particles.n_liquid_particles;
    renderer.render_info.n_boundary_particles = bsph.particles.n_boundary_particles;
    renderer.render_info.n_mobile_particles = bsph.particles.n_mobile_particles;
    renderer.render_info.n_total_particles = bsph.particles.n_total_particles;

    renderer.render_info.x = bsph.particles.x;
    renderer.render_info.y = bsph.particles.y;
    renderer.render_info.z = bsph.particles.z;

    renderer.render_info.vx = bsph.particles.vx;
    renderer.render_info.vy = bsph.particles.vy;
    renderer.render_info.vz = bsph.particles.vz;

    renderer.render_info.Fx = bsph.particles.Fx;
    renderer.render_info.Fy = bsph.particles.Fy;
    renderer.render_info.Fz = bsph.particles.Fz;

    renderer.render_info.rho = bsph.particles.rho;

    renderer.render_info.p = bsph.particles.p;

    renderer.set_grid( &bsph.uniform_grid );

    renderer.init(argc,argv);
    unsigned int it=0;
    bool render = true;

    //load_model_data( 0.15 );

    while(!glfwWindowShouldClose(renderer.getWindow())){
        std::cout << "[" << it * simState.dt << "] sec\n";
        //std::cout<<"iteration "<<it<<std::endl;
        //usleep(10000);
        if(render and it >= startIt && it <= endIt){
            maya.writeToMaya(it,bsph.particles.x,
                                bsph.particles.y,
                                bsph.particles.z,
                                bsph.particles.n_liquid_particles);
            maya.writeBoatInfos(bsph.mobile_mass_center_x,
                                bsph.mobile_mass_center_y,
                                bsph.mobile_mass_center_z,
                                bsph.mobile_angle_phi,
                                bsph.mobile_angle_theta,
                                bsph.mobile_angle_psi,
                                it);
        } 
        bsph.run_step( simState.dt );
        renderer.draw();
        it++;

        //std::cout<<bsph.mobile_mass_center_x<<","<<bsph.mobile_mass_center_y<<","<<bsph.mobile_mass_center_z<<std::endl;
        //std::cout<<bsph.mobile_scale<<std::endl;

    }
#endif


    return 0;
}

void initRendererWithSimInfo(const PCI_SPH &pciSph, OpenGL_Renderer &renderer) {
    renderer.render_info.n_liquid_particles     = pciSph.particles.n_liquid_particles;
    renderer.render_info.n_boundary_particles   = pciSph.particles.n_boundary_particles;
    renderer.render_info.n_mobile_particles     = pciSph.particles.n_mobile_particles;
    renderer.render_info.n_total_particles      = pciSph.particles.n_total_particles;

    renderer.render_info.x = pciSph.particles.x;
    renderer.render_info.y = pciSph.particles.y;
    renderer.render_info.z = pciSph.particles.z;

    renderer.render_info.vx = pciSph.particles.vx;
    renderer.render_info.vy = pciSph.particles.vy;
    renderer.render_info.vz = pciSph.particles.vz;

    renderer.render_info.Fx = pciSph.particles.Fx;
    renderer.render_info.Fy = pciSph.particles.Fy;
    renderer.render_info.Fz = pciSph.particles.Fz;

    renderer.render_info.rho = pciSph.particles.rho;

    renderer.render_info.p = pciSph.particles.p;

    renderer.set_grid( pciSph.uniform_grid );
}







