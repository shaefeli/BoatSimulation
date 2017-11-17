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
    GLuint box_points_VBO, box_indices_VBO;

	GLuint MatrixID;//uniform
	glm::mat4 Projection;
	glm::mat4 View;
	glm::mat4 Model;
	glm::mat4 MVP;

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
        void draw_particles(size_t n_particles, const float *x, const float *y, const float *z);

        void reshape(int x, int y);

};
