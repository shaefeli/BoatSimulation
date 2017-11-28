#include "Basic_SPH_System.h"
#include "OpenGL_Renderer.h"
#include "unistd.h"
#include <iostream>
#define nr_particles 2


int main(int argc, char** argv){

    
    SimState state;
    state.dt = 1e-3;
    state.g  = 9.8;
    state.h  = 0.0625;
    state.k  = 3.5;
    state.mu = 3.5;
    state.rho0 = 1000;
    
    Basic_SPH_System bsph(5000,
                          0., 0., 0.,
                          0.5, 0.5, 0.5,
                          state.h , state.h , state.h );
    bsph.setSimState(state);
    bsph.finilizeInit();

    OpenGL_Renderer renderer;

    renderer.render_info.n_particles = bsph.particles.n_particles;
    renderer.render_info.x = bsph.particles.x;
    renderer.render_info.y = bsph.particles.y;
    renderer.render_info.z = bsph.particles.z;

    renderer.set_grid( &bsph.uniform_grid );

    renderer.init(argc,argv);
    unsigned int it = 0;
    
    while(!glfwWindowShouldClose(renderer.getWindow())){
        std::cout << "[" << it * state.dt << "] sec\n";
        //std::cout<<"iteration "<<it<<std::endl;
//        usleep(10000);
        bsph.run_step( state.dt );
        renderer.draw();
        it++;
    }

    return 0;
}
