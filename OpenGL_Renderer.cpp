#include "OpenGL_Renderer.h"
#include <iostream>
#include <cstdio>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>

#include "./viridis.h"


void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    Render_mode *rm = (Render_mode *)glfwGetWindowUserPointer(window);
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
    {
        glfwSetWindowShouldClose(window, GL_TRUE);
    } 
    else if ( rm != NULL && key == GLFW_KEY_1 && action == GLFW_PRESS)
    {
        *rm = NONE;
    }
    else if ( rm != NULL && key == GLFW_KEY_2 && action == GLFW_PRESS)
    {
        *rm = NEIGHBORS;
    }
    else if ( rm != NULL && key == GLFW_KEY_3 && action == GLFW_PRESS)
    {
        *rm = GRID;
    }
    else if ( rm != NULL && key == GLFW_KEY_4 && action == GLFW_PRESS)
    {
        *rm = FORCE;
    }
    else if ( rm != NULL && key == GLFW_KEY_5 && action == GLFW_PRESS)
    {
        *rm = DENSITY;
    }
    else if ( rm != NULL && key == GLFW_KEY_6 && action == GLFW_PRESS)
    {
        *rm = PRESSURE;
    }
}

GLuint LoadShaders(const char * vertex_file_path,const char * fragment_file_path) {

	// Create the shaders
	GLuint VertexShaderID = glCreateShader(GL_VERTEX_SHADER);
	GLuint FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);

	// Read the Vertex Shader code from the file
	std::string VertexShaderCode;
	std::ifstream VertexShaderStream(vertex_file_path, std::ios::in);
    printf("%s\n",glGetString(GL_VERSION));
	if(VertexShaderStream.is_open()){
		std::string Line = "";
		while(getline(VertexShaderStream, Line))
			VertexShaderCode += "\n" + Line;
		VertexShaderStream.close();
	}else{
		printf("Impossible to open %s. Are you in the right directory ? Don't forget to read the FAQ !\n", vertex_file_path);
		getchar();
		return 0;
	}

	// Read the Fragment Shader code from the file
	std::string FragmentShaderCode;
	std::ifstream FragmentShaderStream(fragment_file_path, std::ios::in);
	if(FragmentShaderStream.is_open()){
		std::string Line = "";
		while(getline(FragmentShaderStream, Line))
			FragmentShaderCode += "\n" + Line;
		FragmentShaderStream.close();
	}

	GLint Result = GL_FALSE;
	int InfoLogLength;


	// Compile Vertex Shader
	printf("Compiling shader : %s\n", vertex_file_path);
	char const * VertexSourcePointer = VertexShaderCode.c_str();
	glShaderSource(VertexShaderID, 1, &VertexSourcePointer , NULL);
	glCompileShader(VertexShaderID);

	// Check Vertex Shader
	glGetShaderiv(VertexShaderID, GL_COMPILE_STATUS, &Result);
	glGetShaderiv(VertexShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if ( InfoLogLength > 0 ){
		std::vector<char> VertexShaderErrorMessage(InfoLogLength+1);
		glGetShaderInfoLog(VertexShaderID, InfoLogLength, NULL, &VertexShaderErrorMessage[0]);
		printf("%s\n", &VertexShaderErrorMessage[0]);
	}



	// Compile Fragment Shader
	printf("Compiling shader : %s\n", fragment_file_path);
	char const * FragmentSourcePointer = FragmentShaderCode.c_str();
	glShaderSource(FragmentShaderID, 1, &FragmentSourcePointer , NULL);
	glCompileShader(FragmentShaderID);

	// Check Fragment Shader
	glGetShaderiv(FragmentShaderID, GL_COMPILE_STATUS, &Result);
	glGetShaderiv(FragmentShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if ( InfoLogLength > 0 ){
		std::vector<char> FragmentShaderErrorMessage(InfoLogLength+1);
		glGetShaderInfoLog(FragmentShaderID, InfoLogLength, NULL, &FragmentShaderErrorMessage[0]);
		printf("%s\n", &FragmentShaderErrorMessage[0]);
	}



	// Link the program
	printf("Linking program\n");
	GLuint ProgramID = glCreateProgram();
	glAttachShader(ProgramID, VertexShaderID);
	glAttachShader(ProgramID, FragmentShaderID);
	glLinkProgram(ProgramID);

	// Check the program
	glGetProgramiv(ProgramID, GL_LINK_STATUS, &Result);
	glGetProgramiv(ProgramID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if ( InfoLogLength > 0 ){
		std::vector<char> ProgramErrorMessage(InfoLogLength+1);
		glGetProgramInfoLog(ProgramID, InfoLogLength, NULL, &ProgramErrorMessage[0]);
		printf("%s\n", &ProgramErrorMessage[0]);
	}

	
	glDetachShader(ProgramID, VertexShaderID);
	glDetachShader(ProgramID, FragmentShaderID);
	
	glDeleteShader(VertexShaderID);
	glDeleteShader(FragmentShaderID);
    printf("correctly created and linked shaders\n");
	return ProgramID;
}

OpenGL_Renderer::OpenGL_Renderer() : ug(NULL)
{
    ;
}

bool OpenGL_Renderer::init( int argc, char** argv)
{
    glutInit(&argc,argv);
    if( !glfwInit() )
    {
            fprintf( stderr, "Failed to initialize GLFW\n" );
                return false;
    }

    //glfwWindowHint(GLFW_SAMPLES, 4);//Antialiasing
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    // Open a window and create its OpenGL context
    //GLFWwindow* window; // this variable is global for simplicity)
    window = glfwCreateWindow( 1024, 768, "Tutorial 01", NULL, NULL);
    if( !window ){
        fprintf( stderr, "Failed to open GLFW window.\n" );
        glfwTerminate();
        exit(EXIT_FAILURE);
    }
    glfwMakeContextCurrent(window); // Initialize GLEW
    glewExperimental=GL_TRUE; // Needed in core profile
    if (glewInit() != GLEW_OK) {
        fprintf(stderr, "Failed to initialize GLEW\n");
        return false;
    }
    
    // Set the keyboard callback so that when we press ESC, it knows what to do.
    glfwSetKeyCallback(window, key_callback);
    glfwSetWindowUserPointer(this->window,(void *)&(this->render_mode));

    programID = LoadShaders("./vertex_shader.vs","./fragment_shader.fs");
    glUseProgram(programID);


	MatrixID = glGetUniformLocation(programID, "MVP");
	color_uniform = glGetUniformLocation(programID, "fragment_color");

	Projection = glm::perspective(glm::radians(45.0f), 4.0f / 3.0f, 0.1f, 100.0f);
	View       = glm::lookAt(
                        glm::vec3(2,1.5,1.5), // Camera is at (4,3,3), in World Space
                        glm::vec3(0,0,0), // and looks at the origin
                        glm::vec3(0,1,0)  // Head is up (set to 0,-1,0 to look upside-down)
                   );
	Model      = glm::mat4(1.0f);
	MVP        = Projection * View * Model;
	//MVP      = glm::mat4(1.0f);

    static const GLfloat g_vertex_buffer_data[] = { 
        //Bottom face
		0.0f, 0.0f, 0.0f,
		1.0f, 0.0f, 0.0f,

		1.0f, 0.0f, 0.0f,
		1.0f, 0.0f, 1.0f,

		1.0f, 0.0f, 1.0f,
		0.0f, 0.0f, 1.0f,

		0.0f, 0.0f, 1.0f,
		0.0f, 0.0f, 0.0f,

        //Upwards face
		0.0f, 1.0f, 0.0f,
		1.0f, 1.0f, 0.0f,

		1.0f, 1.0f, 0.0f,
		1.0f, 1.0f, 1.0f,

		1.0f, 1.0f, 1.0f,
		0.0f, 1.0f, 1.0f,

		0.0f, 1.0f, 1.0f,
		0.0f, 1.0f, 0.0f,

        //Merge them
		0.0f, 0.0f, 0.0f,
		0.0f, 1.0f, 0.0f,

		1.0f, 0.0f, 0.0f,
		1.0f, 1.0f, 0.0f,

		1.0f, 0.0f, 1.0f,
		1.0f, 1.0f, 1.0f,

		0.0f, 0.0f, 1.0f,
		0.0f, 1.0f, 1.0f,


	};
    static const GLfloat g_vertex_buffer_color[] = { 
        //Bottom face
		1.0f, 0.0f, 0.0f,
		1.0f, 0.0f, 0.0f,

		1.0f, 0.0f, 0.0f,
		1.0f, 0.0f, 0.0f,

		1.0f, 0.0f, 0.0f,
		1.0f, 0.0f, 0.0f,

		1.0f, 0.0f, 0.0f,
		1.0f, 0.0f, 0.0f,


		1.0f, 0.0f, 0.0f,
		1.0f, 0.0f, 0.0f,

		1.0f, 0.0f, 0.0f,
		1.0f, 0.0f, 0.0f,

		1.0f, 0.0f, 0.0f,
		1.0f, 0.0f, 0.0f,

		1.0f, 0.0f, 0.0f,
		1.0f, 0.0f, 0.0f,


		1.0f, 0.0f, 0.0f,
		1.0f, 0.0f, 0.0f,

		1.0f, 0.0f, 0.0f,
		1.0f, 0.0f, 0.0f,

		1.0f, 0.0f, 0.0f,
		1.0f, 0.0f, 0.0f,

		1.0f, 0.0f, 0.0f,
		1.0f, 0.0f, 0.0f,
	};


    glPointSize(4.);
	
    
    //Enable depth test
    glEnable(GL_DEPTH_TEST);
    //Accept fragment if it closer to the camera than the former one
    glDepthFunc(GL_LESS); 

    //GLuint box_VAO;
	glGenVertexArrays(1, &box_VAO);
	glBindVertexArray(box_VAO);

    
	//GLuint box_points_VBO;
	glGenBuffers(1, &box_points_VBO);
	glBindBuffer(GL_ARRAY_BUFFER, box_points_VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data), g_vertex_buffer_data, GL_STATIC_DRAW);

	glGenBuffers(1, &box_colors_VBO);
	glBindBuffer(GL_ARRAY_BUFFER, box_colors_VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_color), g_vertex_buffer_color, GL_STATIC_DRAW);
    glBindVertexArray(0);

	glGenVertexArrays(1, &particles_VAO);
	glBindVertexArray(particles_VAO);

	glGenBuffers(1, &particles_VBO);
	glBindBuffer(GL_ARRAY_BUFFER, particles_VBO);

    
    size_t nparticles = render_info.n_particles;
    printf("n_particles: %lu\n",nparticles);
    //nparticles = 5;

    //float *interleaved_data = (float *)malloc(3*nparticles*sizeof(float));
    interleaved_buffer = (float *)malloc(3*nparticles*sizeof(float));
    for( int i = 0; i < nparticles; i++ ) {
        interleaved_buffer[i*3 + 0] = render_info.x[i];
        interleaved_buffer[i*3 + 1] = render_info.y[i];
        interleaved_buffer[i*3 + 2] = render_info.z[i];
    }

	glBufferData(GL_ARRAY_BUFFER, nparticles*3*sizeof(GLfloat), interleaved_buffer, GL_STATIC_DRAW);

	glGenBuffers(1, &colors_VBO);
	glBindBuffer(GL_ARRAY_BUFFER, colors_VBO);

    color_buffer = (float *)malloc(3*nparticles*sizeof(float));
    for( size_t i = 0; i < nparticles; i++ ) {
        color_buffer[3*i]     = 0.;
        color_buffer[3*i+1]   = 0.;
        color_buffer[3*i+2]   = 1.;
    }

    glBufferData(GL_ARRAY_BUFFER, nparticles*3*sizeof(GLfloat), color_buffer, GL_STATIC_DRAW);

    
    glBindVertexArray(0);
    return true;
}



//We don't wanna modify
void OpenGL_Renderer::draw()
{

    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    glUseProgram(programID);
    glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);
    
    draw_particles();
    draw_box();

	glfwSwapBuffers(window);
    glfwPollEvents();
    std::cout<<render_mode<<std::endl;
}

void OpenGL_Renderer::draw_box()
{
        glUniform3f(color_uniform, 1., 0., 0.);
		// Clear the screen
		glBindVertexArray(box_VAO);
        //glEnableVertexArrayAttrib(box_VAO,0);
        glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, box_points_VBO);
		glVertexAttribPointer(
			0,                  // attribute
			3,                  // size
			GL_FLOAT,           // type
			GL_FALSE,           // normalized?
			0,                  // stride
			(void*)0            // array buffer offset
		);

        glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, box_colors_VBO);
		glVertexAttribPointer(
			1,                  // attribute
			3,                  // size
			GL_FLOAT,           // type
			GL_FALSE,           // normalized?
			0,                  // stride
			(void*)0            // array buffer offset
		);




		// Draw the triangle !
		glDrawArrays(GL_LINES, 0, 2*12); // 3 indices starting at 0 -> 1 triangle
        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(1);
		glBindVertexArray(0);
}

void OpenGL_Renderer::draw_particles( )
{
    glUniform3f(color_uniform, 0., 0., 1.);

    //Draw the points
    //glEnableVertexArrayAttrib(particles_VAO,0);
    glBindVertexArray(particles_VAO);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, particles_VBO);

    size_t nparticles = render_info.n_particles;
    //float *interleaved_data = (float *)malloc(3*nparticles*sizeof(float));
    for( int i = 0; i < nparticles; i++ ) {
        interleaved_buffer[i*3 + 0] = render_info.x[i];
        interleaved_buffer[i*3 + 1] = render_info.y[i];
        interleaved_buffer[i*3 + 2] = render_info.z[i];
    }
	glBufferData(GL_ARRAY_BUFFER, nparticles*3*sizeof(GLfloat), interleaved_buffer, GL_STATIC_DRAW);
    //free(interleaved_data);
    
    glVertexAttribPointer(
        0,
        3,
        GL_FLOAT,
        GL_FALSE,
        0,
        (void*)0
    );
    //glfwSetWindowUserPointer
    
    //glEnableVertexArrayAttrib(particles_VAO,1);
    glEnableVertexAttribArray(1);
	glBindBuffer(GL_ARRAY_BUFFER, colors_VBO);


    //Set the color depending on the render mode
    switch(this->render_mode) {
        case NEIGHBORS:
            set_neighbor_color(0);
            break;
        case GRID:
            set_grid_color();
            break;
        case SPEED:
            break;
        case FORCE:
            break;
        case DENSITY:
            set_density_color();
            break;
        case PRESSURE:
            set_pressure_color();
            break;
        default:

            break;
    }
    //set_grid_color();
    //set_neighbor_color(0);

    glBufferData(GL_ARRAY_BUFFER, nparticles*3*sizeof(GLfloat), color_buffer, GL_STATIC_DRAW);
 
    glVertexAttribPointer(
        1,
        3,
        GL_FLOAT,
        GL_FALSE,
        0,
        (void*)0
    );
    

    glDrawArrays(GL_POINTS, 0, nparticles*3);

    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);


}


void OpenGL_Renderer::reshape(int x, int y)
{
    if (y == 0 || x == 0) return;  //Nothing is visible then, so return
    //Set a new projection matrix
    //Angle of view:40 degrees
    //Near clipping plane distance: 0.5
    //Far clipping plane distance: 20.0
     
}




void OpenGL_Renderer::set_grid_color()
{
    size_t nparticles = render_info.n_particles;
    //Uniform_Grid ug(0,0,0, 1,1,1, 0.1,0.1,0.1);
    for( size_t i = 0; i < nparticles; i++ ) {
        size_t gi,gj,gk;
        ug->query_cell(render_info.x[i],render_info.y[i],render_info.z[i],gi,gj,gk);
        color_buffer[3*i]   = float(gi)/ug->n_cells_x;
        color_buffer[3*i+1] = float(gj)/ug->n_cells_y;
        color_buffer[3*i+2] = float(gk)/ug->n_cells_z;
    }
}


void OpenGL_Renderer::set_grid( Uniform_Grid *ug )
{
    this->ug = ug;
}

void OpenGL_Renderer::set_neighbor_color( size_t particle_index )
{
    //Uniform_Grid ug(0,0,0, 1,1,1, 0.1,0.1,0.1);
    //ug.build(render_info.x,render_info.y,render_info.z,render_info.n_particles);
   

    for( size_t i = 0; i < render_info.n_particles; i++ ) {
        color_buffer[3*i]   = 0.;
        color_buffer[3*i+1] = 0.;
        color_buffer[3*i+2] = 1.;
    }

    if( ug == NULL ) return;

    std::vector<size_t> neighbor_cells;
    ug->query_neighbors(render_info.x[particle_index],
                       render_info.y[particle_index],
                       render_info.z[particle_index],
                       neighbor_cells);

    int alternate = 0;
    for( size_t i = 0; i < neighbor_cells.size(); i++ ) {
        assert(neighbor_cells[i] < ug->n_cells);
        size_t *cell = ug->cells[neighbor_cells[i]];
        size_t cell_size = ug->cell_size[neighbor_cells[i]];
        for( int j = 0; j < cell_size; j++ ) {
            size_t particle_pos = cell[j];
            if( alternate == 0 ) {
                color_buffer[3*particle_pos]   = 1.;
                color_buffer[3*particle_pos+1] = 0.;
                color_buffer[3*particle_pos+2] = 0.;
            } else if( alternate == 1 ){
                color_buffer[3*particle_pos]   = 0.;
                color_buffer[3*particle_pos+1] = 1.;
                color_buffer[3*particle_pos+2] = 0.;
            } else {
                color_buffer[3*particle_pos]   = 1.;
                color_buffer[3*particle_pos+1] = 0.;
                color_buffer[3*particle_pos+2] = 1.;
            }
        }
        alternate = (alternate+1)%2;
    }
    color_buffer[3*particle_index]   = 1.;
    color_buffer[3*particle_index+1] = 1.;
    color_buffer[3*particle_index+2] = 0.;

    ug->clean_up();

}



void OpenGL_Renderer::set_density_color()
{
    double min = 2000.;
    double max = 3500.;
    size_t nparticles = render_info.n_particles;
    for( size_t i = 0; i < nparticles; i++ ) {
        //std::cout<<(render_info.rho[i]-min)/(max-min)<<std::endl;
        //std::cout<<render_info.rho[i]<<std::endl;
        float val = std::max(0.,std::min(1.,(render_info.rho[i]-min)/(max-min)));
        int index = val*256.;
        //std::cout<<val<<std::endl;
        //std::cout<<index<<std::endl;
        color_buffer[3*i]   = viridis_data[index][0];
        color_buffer[3*i+1] = viridis_data[index][1];
        color_buffer[3*i+2] = viridis_data[index][2];
    }
}


void OpenGL_Renderer::set_pressure_color()
{
    double min = 5000.;
    double max = 15000.;
    size_t nparticles = render_info.n_particles;
    for( size_t i = 0; i < nparticles; i++ ) {
        //std::cout<<render_info.p[i]<<std::endl;
        //std::cout<<(render_info.p[i]-min)/(max-min)<<std::endl;
        //std::cout<<render_info.p[i]<<std::endl;
        float val = std::max(0.,std::min(1.,(render_info.p[i]-min)/(max-min)));
        int index = val*256.;
        //std::cout<<val<<std::endl;
        //std::cout<<index<<std::endl;
        color_buffer[3*i]   = viridis_data[index][0];
        color_buffer[3*i+1] = viridis_data[index][1];
        color_buffer[3*i+2] = viridis_data[index][2];
    }
}







