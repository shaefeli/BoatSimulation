#include "OpenGL_Renderer.h"
#include "Basic_SPH_System.h"
#include "Maya_Interface.h"
#include "unistd.h"
#include "Particle_Generator.h"
#include <iostream>


int main(int argc, char** argv){

    int endIt = 400;
    int startIt = 301;

    SimState state;
    state.dt = 5e-4;
    //state.dt = 1e-5;
    state.g  = 9.8;
    state.h  = 0.0625;
    state.k  = 3.5;
    state.mu = 3.5;
    state.rho0 = 1000;
    
    Basic_SPH_System bsph(
                          0., 0., 0.,
                          1., 1., 1.,
                          state.h , state.h , state.h );
    bsph.setSimState(state);
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
        std::cout << "[" << it * state.dt << "] sec\n";
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
        bsph.run_step( state.dt );
        renderer.draw();
        it++;
        //std::cout<<bsph.mobile_mass_center_x<<","<<bsph.mobile_mass_center_y<<","<<bsph.mobile_mass_center_z<<std::endl;
        //std::cout<<bsph.mobile_scale<<std::endl;

    }
    std::cout<<"closing"<<std::endl;

    return 0;
}







