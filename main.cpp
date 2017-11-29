#include "OpenGL_Renderer.h"
#include "Basic_SPH_System.h"
#include "Maya_Interface.h"
#include "unistd.h"
#include <iostream>


int main(int argc, char** argv){

    
    SimState state;
    state.dt = 5e-4;
    state.g  = 9.8;
    state.h  = 0.0625;
    state.k  = 3.5;
    state.mu = 3.5;
    state.rho0 = 1000;
    
    Basic_SPH_System bsph(5000,
                          0., 0., 0.,
                          0.5, .5, .5,
                          state.h , state.h , state.h );
    bsph.setSimState(state);
    bsph.finilizeInit();

    OpenGL_Renderer renderer;
    Maya_Interface maya(5000);

    renderer.render_info.n_particles = bsph.particles.n_particles;
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

    while(!glfwWindowShouldClose(renderer.getWindow())){
        std::cout << "[" << it * state.dt << "] sec\n";
        //std::cout<<"iteration "<<it<<std::endl;
//        usleep(10000);
        if(it > 300 && it < 400){
            maya.writeToMaya(it,bsph.particles.x,bsph.particles.y,bsph.particles.z, bsph.particles.n_particles);
        } 
        bsph.run_step( state.dt );
        renderer.draw();
        it++;
    }

    return 0;
}
