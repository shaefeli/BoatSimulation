#include <GL/glew.h>
#include <GL/glut.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <cstddef>


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

    float *interleaved_buffer;
    float *color_buffer;

    public:
 
        typedef struct {
                            size_t n_particles;
                            const float *x;
                            const float *y;
                            const float *z;
        } rendering_info_t;
       
        rendering_info_t render_info;

        bool init( int argc, char** argv);
        void draw();

        void draw_box();
        void draw_particles();
    
        //Function to set the color according to parameter
        void set_grid_color();
        void set_neighbor_color( size_t particle_index );

        void reshape(int x, int y);

};
