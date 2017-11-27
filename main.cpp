#include "OpenGL_Renderer.h"
#include "Basic_SPH_System.h"
#include "Maya_Interface.h"
#include "unistd.h"
#include <iostream>


int main(int argc, char** argv){

    size_t nr_particles =1000;
    Basic_SPH_System bsph(nr_particles);
    OpenGL_Renderer renderer;
    Maya_Interface maya(nr_particles);

    renderer.render_info.n_particles = bsph.particles.n_particles;
    renderer.render_info.x = bsph.particles.x;
    renderer.render_info.y = bsph.particles.y;
    renderer.render_info.z = bsph.particles.z;

    renderer.init(argc,argv);
    unsigned int it=0;
    while(true) {
        usleep(10000);
        bsph.run_step( 0.01 );
        //This if statement is temporary so it doesn't write to file for every iteration
        if(it == 0 || it==1 || it==2){
            maya.writeToMaya(it,bsph.particles.x,bsph.particles.y,bsph.particles.z, bsph.particles.n_particles);
        }    
        renderer.draw();
        it++;

    }

    return 0;
}
