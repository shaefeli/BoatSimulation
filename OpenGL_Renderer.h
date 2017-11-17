#include <GL/glew.h>
#include <GL/glut.h>
#include <GLFW/glfw3.h>
#include <cstddef>

class OpenGL_Renderer
{
    private:
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
