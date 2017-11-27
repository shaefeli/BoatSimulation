#include "Basic_SPH_System.h"
#include "OpenGL_Renderer.h"
#include "unistd.h"
#include <iostream>
#define nr_particles 2


int main(int argc, char** argv){

    Basic_SPH_System bsph(1000, 0., 0., 0., 1., 1., 1., 0.1, 0.1, 0.1);
    SimState state;
    state.dt = 1e-4;
    state.g  = 9.8;
    state.h  = 5e-2;
    state.k  = 1e3;
    state.mu = 0.1;
    state.rho0 = 1000;
    bsph.setSimState(state);
    bsph.finilizeInit();

    OpenGL_Renderer renderer;

    renderer.render_info.n_particles = bsph.particles.n_particles;
    renderer.render_info.x = bsph.particles.x;
    renderer.render_info.y = bsph.particles.y;
    renderer.render_info.z = bsph.particles.z;

    renderer.set_grid( &bsph.uniform_grid );

    renderer.init(argc,argv);
    unsigned int it;
    while(!glfwWindowShouldClose(renderer.getWindow())){
        //std::cout<<"iteration "<<it<<std::endl;
//        usleep(10000);
        bsph.run_step( state.dt );
        renderer.draw();
        it++;
    }

    return 0;
}
