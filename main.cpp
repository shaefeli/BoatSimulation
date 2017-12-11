#include "OpenGL_Renderer.h"
#include "Basic_SPH_System.h"
#include "PCI_SPH.h"
#include "SupportingStructures.h"
#include "Maya_Interface.h"
#include "unistd.h"
#include "Particle_Generator.h"
#include <iostream>

#define METHOD_BASIC_SPH 0
#define METHOD_PCI_SPH   1

void initRendererWithSimInfo(const PCI_SPH &pciSph, OpenGL_Renderer &renderer);

int main(int argc, char** argv){

    
    SimState simState;
    simState.dt = 5e-4;
    //simState.dt = 1e-5;
    simState.g  = 9.8;
    simState.h  = 0.0625;
    simState.k  = 3.5;
    simState.mu = 3.5;
    simState.rho0 = 1000;

#if METHOD_PCI_SPH


    BoundaryBox bBox{
            .x1 = 0.0f, .x2 = 1.0f,
            .y1 = 0.0f, .y2 = 1.0f,
            .z1 = 0.0f, .z2 = 1.0f
    };

    ParticlesInitialSpawningBox iBox{
            .x1 = 0.25f,  .x2 = .75f,
            .y1 = 0.001f, .y2 = .5f,
            .z1 = 0.25f,  .z2 = .75f,
            .spawningRadius = simState.h * 0.7f
    };

    UniformGridSplit gridSplit{
            .cells_x = simState.h,
            .cells_y = simState.h,
            .cells_z = simState.h
    };

    PCI_SPH pciSph(bBox,iBox,simState,gridSplit);


    OpenGL_Renderer renderer;
    initRendererWithSimInfo(pciSph, renderer);
    renderer.init(argc,argv);

    while(!glfwWindowShouldClose(renderer.getWindow())){
        std::cout << "[" << pciSph.getCurrentTime() << "] sec\n";
        pciSph.run_step( );
        renderer.draw();
    }


#endif


#if METHOD_BASIC_SPH
    Basic_SPH_System bsph(
                          0., 0., 0.,
                          1., 1., 1.,
                          simState.h , simState.h , simState.h );
    bsph.setSimState(simState);
    bsph.finilizeInit();

    OpenGL_Renderer renderer;
    Maya_Interface maya(5000);

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
    bool render = false;

    //load_model_data( 0.15 );

    while(!glfwWindowShouldClose(renderer.getWindow())){
        std::cout << "[" << it * simState.dt << "] sec\n";
        //std::cout<<"iteration "<<it<<std::endl;
        //usleep(10000);
//        if(render and it > 300 && it < 400){
//            maya.writeToMaya(it,bsph.particles.x,
//                                bsph.particles.y,
//                                bsph.particles.z,
//                                bsph.particles.n_liquid_particles);
//        }
        bsph.run_step( simState.dt );
        renderer.draw();
        it++;
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

    renderer.set_grid( pciSph.uniform_grid.get() );
}







