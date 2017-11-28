#ifndef _PBS_OPENGL_RENDERER_
#define _PBS_OPENGL_RENDERER_
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <GL/glut.h>
#include <cstddef>
#include "Uniform_Grid.h"


//enum Render_Mode = { NEIGHBORS, GRID, SPEED, FORCE, DENSITY, NONE };


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

    //Render_mode render_mode = NONE;

    public:
    
        typedef struct {
                            size_t n_particles;
                            const float *x;
                            const float *y;
                            const float *z;
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

        void reshape(int x, int y);
        
        void set_grid( Uniform_Grid *ug );
};

#endif
