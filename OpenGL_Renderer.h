#ifndef _PBS_OPENGL_RENDERER_
#define _PBS_OPENGL_RENDERER_
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/ext.hpp>

#ifdef __APPLE__
    #include <GLUT/glut.h>
#else
    #include <GL/glut.h>
#endif

#include <cstddef>
#include "Uniform_Grid.h"



typedef enum { NONE, NEIGHBORS, GRID, SPEED, FORCE, DENSITY, PRESSURE } Render_mode;


class OpenGL_Renderer
{
    private:
    GLFWwindow* window; // this variable is global for simplicity)
    GLuint VertexArrayID;
    GLuint vertexbuffer;
    GLuint programID;

    GLuint box_VAO;
    GLuint box_points_VBO,box_colors_VBO;

    GLuint particles_VAO;
    GLuint particles_VBO, colors_VBO;

    GLuint element_index_VBO, element_positions_VBO;

    //uniforms
    GLuint MatrixID;
    GLuint color_uniform;

    glm::mat4 Projection;
    glm::mat4 View;
    glm::mat4 Model;
    glm::mat4 MVP;


    Uniform_Grid *ug;

    float *interleaved_buffer;
    float *color_buffer;


    public:
    bool display_boundary = false;
    Render_mode render_mode = NONE;
    
        typedef struct {
                            size_t n_liquid_particles;
                            size_t n_total_particles;

                            const float *x;
                            const float *y;
                            const float *z;

                            const float *vx;
                            const float *vy;
                            const float *vz;

                            const float *Fx;
                            const float *Fy;
                            const float *Fz;
                            
                            const float *rho;   //density
                            const float *p;     //pressure

        } rendering_info_t;

        OpenGL_Renderer();
        rendering_info_t render_info;

    
        GLFWwindow* getWindow(){ return this->window;};
        bool init( int argc, char** argv);
        void draw();

        void draw_box();
        void draw_particles();
    
        //Function to set the color according to parameter
        void set_grid_color();
        void set_neighbor_color( size_t particle_index );
        void set_density_color();
        void set_pressure_color();

        void reshape(int x, int y);
        
        void set_grid( Uniform_Grid *ug );
};

#endif
