#include "OpenGL_Renderer.h"
#include "Basic_SPH_System.h"
#include "unistd.h"
#include <iostream>
#define nr_particles 2


int main(int argc, char** argv){
    OpenGL_Renderer renderer;
    renderer.init(argc,argv);

    Basic_SPH_System bsph(30);

    renderer.render_info.n_particles = bsph.particles.n_particles;
    renderer.render_info.x = bsph.particles.x;
    renderer.render_info.y = bsph.particles.y;
    renderer.render_info.z = bsph.particles.z;
    unsigned int it;
    while(true) {
        //std::cout<<"iteration "<<it<<std::endl;
        usleep(100000);
        bsph.run_step( 0.01 );
        renderer.draw();
        it++;

    }

    return 0;
}
