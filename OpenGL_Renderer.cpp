#include "OpenGL_Renderer.h"
#include <iostream>
#include <cstdio>
#include <fstream>
#include <limits>
#include <cstdlib>
#include <string>
#include <vector>
#include <glm/gtx/transform.hpp>

#include "./viridis.h"


void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    OpenGL_Renderer *oglr = (OpenGL_Renderer *)glfwGetWindowUserPointer(window);
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
    {
        glfwSetWindowShouldClose(window, GL_TRUE);
    } 
    else if ( oglr != NULL && key == GLFW_KEY_1 && action == GLFW_PRESS)
    {
        oglr->render_mode = NONE;
    }
    else if ( oglr != NULL && key == GLFW_KEY_2 && action == GLFW_PRESS)
    {
        oglr->render_mode = NEIGHBORS;
    }
    else if ( oglr != NULL && key == GLFW_KEY_3 && action == GLFW_PRESS)
    {
        oglr->render_mode = GRID;
    }
    else if ( oglr != NULL && key == GLFW_KEY_4 && action == GLFW_PRESS)
    {
        oglr->render_mode = FORCE;
    }
    else if ( oglr != NULL && key == GLFW_KEY_5 && action == GLFW_PRESS)
    {
        oglr->render_mode = DENSITY;
    }
    else if ( oglr != NULL && key == GLFW_KEY_6 && action == GLFW_PRESS)
    {
        oglr->render_mode = PRESSURE;
    }
    else if ( oglr != NULL && key == GLFW_KEY_0 && action == GLFW_PRESS)
    {
        std::cout<<"display boudns switched"<<std::endl;
        oglr->display_boundary = not oglr->display_boundary;
    }
    else if ( oglr != NULL && key == GLFW_KEY_9 && action == GLFW_PRESS)
    {
        std::cout<<"display model switched"<<std::endl;
        oglr->display_mobile = not oglr->display_mobile;
    }
    else if ( oglr != NULL && key == GLFW_KEY_A && action == GLFW_PRESS)
    {
        oglr->xcam -= 0.25;
        oglr->View = glm::lookAt(
          glm::vec3(oglr->xcam,oglr->ycam,oglr->zcam), // Camera is at (4,3,3), in World Space
          glm::vec3(0,0,0), // and looks at the origin
          glm::vec3(0,1,0)  // Head is up (set to 0,-1,0 to look upside-down)
               );
	    oglr->MVP = oglr->Projection * oglr->View * oglr->Model;
    }
    else if ( oglr != NULL && key == GLFW_KEY_D && action == GLFW_PRESS)
    {
        oglr->xcam += 0.25;
        oglr->View = glm::lookAt(
          glm::vec3(oglr->xcam,oglr->ycam,oglr->zcam), // Camera is at (4,3,3), in World Space
          glm::vec3(0,0,0), // and looks at the origin
          glm::vec3(0,1,0)  // Head is up (set to 0,-1,0 to look upside-down)
               );
	    oglr->MVP = oglr->Projection * oglr->View * oglr->Model;
    }
    else if ( oglr != NULL && key == GLFW_KEY_W && action == GLFW_PRESS)
    {
        oglr->zcam += 0.25;
        oglr->View = glm::lookAt(
          glm::vec3(oglr->xcam,oglr->ycam,oglr->zcam), // Camera is at (4,3,3), in World Space
          glm::vec3(0,0,0), // and looks at the origin
          glm::vec3(0,1,0)  // Head is up (set to 0,-1,0 to look upside-down)
               );
	    oglr->MVP = oglr->Projection * oglr->View * oglr->Model;
    }
    else if ( oglr != NULL && key == GLFW_KEY_S && action == GLFW_PRESS)
    {
        oglr->zcam -= 0.25;
        oglr->View = glm::lookAt(
          glm::vec3(oglr->xcam,oglr->ycam,oglr->zcam), // Camera is at (4,3,3), in World Space
          glm::vec3(0,0,0), // and looks at the origin
          glm::vec3(0,1,0)  // Head is up (set to 0,-1,0 to look upside-down)
               );
	    oglr->MVP = oglr->Projection * oglr->View * oglr->Model;
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
    //glfwSetWindowUserPointer(this->window,(void *)&(this->render_mode));
    glfwSetWindowUserPointer(this->window,(void *)this);

    programID = LoadShaders("./vertex_shader.vs","./fragment_shader.fs");
    glUseProgram(programID);


	MatrixID = glGetUniformLocation(programID, "MVP");
	color_uniform = glGetUniformLocation(programID, "fragment_color");

	Projection = glm::perspective(glm::radians(45.0f), 4.0f / 3.0f, 0.1f, 100.0f);
	View       = glm::lookAt(
                        glm::vec3(xcam,ycam,zcam), // Camera is at (4,3,3), in World Space
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

    
    size_t nparticles = render_info.n_total_particles;
    printf("n_liquid_particles: %lu\n",nparticles);
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
    
    //Element to draw
	glGenBuffers(1, &element_index_VBO);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, element_positions_VBO);


	glGenBuffers(1, &element_normals_VBO);
	glBindBuffer(GL_ARRAY_BUFFER, element_normals_VBO);


    static const GLfloat minimal_element_vertices[] = { 
        //Face 1
        -1.0f,-1.0f, 1.0f, // 0
         1.0f,-1.0f, 1.0f, // 1
        -1.0f, 1.0f, 1.0f, // 2
         1.0f, 1.0f, 1.0f, // 3

        -1.0f,-1.0f,-1.0f, // 4
         1.0f,-1.0f,-1.0f, // 5
        -1.0f, 1.0f,-1.0f, // 6
         1.0f, 1.0f,-1.0f, // 7


        -1.0f, 1.0f,-1.0f, // 8
        -1.0f, 1.0f, 1.0f, // 9
         1.0f, 1.0f,-1.0f, // 10
         1.0f, 1.0f, 1.0f, // 11

        -1.0f,-1.0f,-1.0f, // 12
        -1.0f,-1.0f, 1.0f, // 13
         1.0f,-1.0f,-1.0f, // 14
         1.0f,-1.0f, 1.0f, // 15


         1.0f,-1.0f,-1.0f, // 16
         1.0f,-1.0f, 1.0f, // 17
         1.0f, 1.0f,-1.0f, // 18
         1.0f, 1.0f, 1.0f, // 19

        -1.0f,-1.0f,-1.0f, // 20
        -1.0f,-1.0f, 1.0f, // 21
        -1.0f, 1.0f,-1.0f, // 22
        -1.0f, 1.0f, 1.0f, // 23

    };

    static const GLfloat minimal_element_normals[] = {
         0.0f, 0.0f, 1.0f,
         0.0f, 0.0f, 1.0f,
         0.0f, 0.0f, 1.0f,
         0.0f, 0.0f, 1.0f,

         0.0f, 0.0f,-1.0f,
         0.0f, 0.0f,-1.0f,
         0.0f, 0.0f,-1.0f,
         0.0f, 0.0f,-1.0f,
        
         0.0f, 1.0f, 0.0f,
         0.0f, 1.0f, 0.0f,
         0.0f, 1.0f, 0.0f,
         0.0f, 1.0f, 0.0f,

         0.0f,-1.0f, 0.0f,
         0.0f,-1.0f, 0.0f,
         0.0f,-1.0f, 0.0f,
         0.0f,-1.0f, 0.0f,

        -1.0f, 0.0f, 0.0f,
        -1.0f, 0.0f, 0.0f,
        -1.0f, 0.0f, 0.0f,
        -1.0f, 0.0f, 0.0f,

         1.0f, 0.0f, 0.0f,
         1.0f, 0.0f, 0.0f,
         1.0f, 0.0f, 0.0f,
         1.0f, 0.0f, 0.0f,
    };

    static const GLuint minimal_element_indices[] = {
        0,1,2,
        1,3,2,
        
        4,5,6,
        5,7,6,
        
        8,9,10,
        9,11,10,
        
        12,13,14,
        13,15,14,
        
        16,17,18,
        17,19,18,

        20,21,22,
        21,23,22,
    };


	glGenBuffers(1, &element_positions_VBO);
	glBindBuffer(GL_ARRAY_BUFFER, element_positions_VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(minimal_element_vertices), minimal_element_vertices, GL_STATIC_DRAW);

	glGenBuffers(1, &element_index_VBO);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, element_index_VBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(minimal_element_indices), minimal_element_indices, GL_STATIC_DRAW);
    //Element drawn

	glGenBuffers(1, &element_normals_VBO);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, element_normals_VBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(minimal_element_normals), minimal_element_normals, GL_STATIC_DRAW);
   



    glBindVertexArray(0);
    return true;
}



//We don't wanna modify
void OpenGL_Renderer::draw()
{

    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glClearColor(0.7,0.7,0.7,0.);

    glUseProgram(programID);
    //////
	//MVP      = glm::mat4(1.0f);
    //7////
    glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);
    
    //draw_particles();
    draw_box();
    //draw_element(0.2,0.3,0.8);
    draw_particles_elements();

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

    size_t nparticles;
    if(display_boundary) nparticles = render_info.n_liquid_particles;
    else nparticles = render_info.n_total_particles;

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
    glDisableVertexAttribArray(2);




}


void OpenGL_Renderer::draw_element( float x, float y, float z )
{
    
    //Transform
    Projection = glm::perspective(glm::radians(45.0f), 4.0f / 3.0f, 0.1f, 100.0f);
	View       = glm::lookAt(
                        glm::vec3(xcam,ycam,zcam), // Camera is at (4,3,3), in World Space
                        glm::vec3(0,0,0), // and looks at the origin
                        glm::vec3(0,1,0)  // Head is up (set to 0,-1,0 to look upside-down)
                   );
	//Model      = glm::mat4(1.0f);
    
    glm::mat4 Rotation       = glm::mat4(0.01f);
    Rotation[3][3] = 1.;
    

    glm::mat4 Translation = glm::translate(glm::vec3(x,y,z));

	Model      = Translation*Rotation*glm::mat4(1.0f);
	MVP        = Projection * View * Model;
    glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);
    




    //Draw element
   
    glBindVertexArray(particles_VAO);
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(2);
    glBindBuffer(GL_ARRAY_BUFFER, element_positions_VBO);

    
    glVertexAttribPointer(
        0,
        3,
        GL_FLOAT,
        GL_FALSE,
        0,
        (void*)0
    );

    glVertexAttribPointer(
        2,
        3,
        GL_FLOAT,
        GL_FALSE,
        0,
        (void*)0
    );


    //glfwSetWindowUserPointer
    

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, element_index_VBO);

    //glDrawArrays(GL_TRIANGLES, 0, 3);
    
    glDrawElements(
            GL_TRIANGLES,
            12*3,
            GL_UNSIGNED_INT,
            (void*)0
            );


    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(2);

	MVP        = Projection * View;
    glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);

}



void OpenGL_Renderer::draw_particles_elements( )
{
    glUniform3f(color_uniform, 0., 0., 1.);

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



    //Draw the points
    //glEnableVertexArrayAttrib(particles_VAO,0);

    size_t part_offset = 0;
    //if(display_boundary) nparticles = render_info.n_liquid_particles;
    //else nparticles = render_info.n_total_particles;

    //float *interleaved_data = (float *)malloc(3*nparticles*sizeof(float));
    
    for( int i = 0; i < render_info.n_liquid_particles; i++ ) {
        glUniform3f(color_uniform, color_buffer[3*i], color_buffer[3*i+1], color_buffer[3*i+2]);
        draw_element(render_info.x[i],render_info.y[i],render_info.z[i]);
    }
    part_offset += render_info.n_liquid_particles;
    //std::cout<<"total parts:"<<render_info.n_total_particles<<std::endl;
    //std::cout<<"liq parts:"<<render_info.n_liquid_particles<<std::endl;
    //std::cout<<"parts:"<<part_offset<<std::endl;
    if( not display_boundary ) {;
        //std::cout<<"boundary parts:"<<render_info.n_boundary_particles<<std::endl;
        for( int i = part_offset; i < part_offset + render_info.n_boundary_particles; i++ ) {
            glUniform3f(color_uniform, color_buffer[3*i], color_buffer[3*i+1], color_buffer[3*i+2]);
            draw_element(render_info.x[i],render_info.y[i],render_info.z[i]);
        }
    }
    part_offset += render_info.n_boundary_particles;
    if( not display_mobile ) {;
        std::cout<<"mobile parts:"<<render_info.n_mobile_particles<<std::endl;
        for( int i = part_offset; i < part_offset + render_info.n_mobile_particles; i++ ) {
            glUniform3f(color_uniform, color_buffer[3*i], color_buffer[3*i+1], color_buffer[3*i+2]);
            //glUniform3f(color_uniform, 1., 0., 0.);
            draw_element(render_info.x[i],render_info.y[i],render_info.z[i]);
            //std::cout<<"mobile parts:"<<render_info.x[i]<<","<<render_info.y[i]<<","<<render_info.z[i]<<std::endl;
        }
    }

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
    size_t nparticles;
    if(display_boundary) nparticles = render_info.n_liquid_particles;
    else nparticles = render_info.n_total_particles;
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
    //ug.build(render_info.x,render_info.y,render_info.z,render_info.n_liquid_particles);
   
    size_t nparticles = render_info.n_total_particles;
    //if(display_boundary) nparticles = render_info.n_liquid_particles;
    //else nparticles = render_info.n_total_particles;

    for( size_t i = 0; i < nparticles; i++ ) {
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
    size_t nparticles = render_info.n_total_particles;
    //if(display_boundary) nparticles = render_info.n_liquid_particles;
    //else nparticles = render_info.n_total_particles;

    //double min = 2000.;
    //double max = 3500.;
    
    float min = std::numeric_limits<float>::max();
    float max = std::numeric_limits<float>::min();
    for( size_t i = 0; i < nparticles; i++ ) {
        min = std::min(min,render_info.rho[i]);
        max = std::max(max,render_info.rho[i]);
    }



    for( size_t i = 0; i < nparticles; i++ ) {
        //std::cout<<(render_info.rho[i]-min)/(max-min)<<std::endl;
        //std::cout<<render_info.rho[i]<<std::endl;
        float val = std::max(float(0.),std::min(float(1.),(render_info.rho[i]-min)/(max-min)));
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
    size_t nparticles = render_info.n_total_particles;
    //if(display_boundary) nparticles = render_info.n_liquid_particles;
    //else nparticles = render_info.n_total_particles;
    
    float min = std::numeric_limits<float>::max();
    float max = std::numeric_limits<float>::min();
    for( size_t i = 0; i < nparticles; i++ ) {
        min = std::min(min,render_info.p[i]);
        max = std::max(max,render_info.p[i]);
    }

    for( size_t i = 0; i < nparticles; i++ ) {
        //std::cout<<render_info.p[i]<<std::endl;
        //std::cout<<(render_info.p[i]-min)/(max-min)<<std::endl;
        //std::cout<<render_info.p[i]<<std::endl;
        float val = std::max(float(0.),std::min(float(1.),(render_info.p[i]-min)/(max-min)));
        int index = val*256.;
        //std::cout<<val<<std::endl;
        //std::cout<<index<<std::endl;
        color_buffer[3*i]   = viridis_data[index][0];
        color_buffer[3*i+1] = viridis_data[index][1];
        color_buffer[3*i+2] = viridis_data[index][2];
    }
}







