#include "OpenGL_Renderer.h"

void OpenGL_Renderer::init( int argc, char** argv)
{

    glutInit(&argc, argv);
    //we initizlilze the glut. functions
    glutInitDisplayMode(GLUT_SINGLE|GLUT_RGB);
    glutInitWindowPosition(100, 100);
    glutCreateWindow(argv[1]);

    glClearColor(0,0,0,0);

    //glutDisplayFunc(draw);
    //glutReshapeFunc(reshape);
    //Set the function for the animation.
    //glutIdleFunc(animation);
    //glutMainLoop();

}



//We don't wanna modify
void OpenGL_Renderer::draw()
{
    draw_particles(render_info.n_particles, render_info.x, render_info.y, render_info.z);
    draw_box();
}

void OpenGL_Renderer::draw_box()
{
    glMatrixMode(GL_MODELVIEW);
    // clear the drawing buffer.
    glClear(GL_COLOR_BUFFER_BIT);
    glLoadIdentity();
    glTranslatef(0.0,0.0,-10.5);
    glRotated(45,1,1,1);

    glLineWidth(2.5);
    glColor3f(1.0, 0.0, 0.0);
    glBegin(GL_LINES);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(2, 0, 0);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0, 2, 0);
    glVertex3f(0.0, 2.0, 0.0);
    glVertex3f(2, 2, 0);
    glVertex3f(2, 2, 0.0);
    glVertex3f(2, 0, 0);
    glVertex3f(2, 0, 0.0);
    glVertex3f(0, 0, 0);
    glEnd();

    glBegin(GL_LINES);
    glVertex3f(0, 0, 0.0);
    glVertex3f(0, 0, 2);
    glEnd();

    glBegin(GL_LINES);
    glVertex3f(2, 0, 0.0);
    glVertex3f(2, 0, 2);
    glEnd();

     glBegin(GL_LINES);
    glVertex3f(0, 2, 0.0);
    glVertex3f(0, 2, 2);
    glEnd();

     glBegin(GL_LINES);
    glVertex3f(2, 2, 0.0);
    glVertex3f(2, 2, 2);
    glEnd();

    glBegin(GL_LINES);
    glVertex3f(0.0, 0.0, 2);
    glVertex3f(2, 0, 2);
    glVertex3f(0.0, 0.0, 2);
    glVertex3f(0, 2, 2);
    glVertex3f(0.0, 2.0, 2);
    glVertex3f(2, 2, 2);
    glVertex3f(2, 2, 2);
    glVertex3f(2, 0, 2);
    glVertex3f(2, 0, 2);
    glVertex3f(0, 0, 2);
    glEnd();
 
}

void OpenGL_Renderer::draw_particles(   size_t n_particles,
                                        const float *x,
                                        const float *y,
                                        const float *z
                                    )
{
    for(int i=0; i<n_particles;i++){
        glColor3f(0.0, 1.0, 0.0);
        glPushMatrix();
        float posX = x[i];
        float posY = y[i];
        float posZ = z[i];
        glTranslatef(posX,posY,posZ);
        glutSolidSphere(0.1,50.0,50.0);
        glPopMatrix();
    }

    glFlush();
}


void OpenGL_Renderer::reshape(int x, int y)
{
    if (y == 0 || x == 0) return;  //Nothing is visible then, so return
    //Set a new projection matrix
    glMatrixMode(GL_PROJECTION);  
    glLoadIdentity();
    //Angle of view:40 degrees
    //Near clipping plane distance: 0.5
    //Far clipping plane distance: 20.0
     
    gluPerspective(40.0,(GLdouble)x/(GLdouble)y,0.5,20.0);
    glMatrixMode(GL_MODELVIEW);
    glViewport(0,0,x,y);  //Use the whole window for rendering
}




void animation(void)
{
    //Main function where we update the values of the particles
     //sph->update();

     //draw new frame
     //DrawSPH();
}










