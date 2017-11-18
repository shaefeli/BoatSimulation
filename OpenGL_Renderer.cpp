#include "OpenGL_Renderer.h"
#include <cstdio>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>



GLuint LoadShaders(const char * vertex_file_path,const char * fragment_file_path) {

	// Create the shaders
	GLuint VertexShaderID = glCreateShader(GL_VERTEX_SHADER);
	GLuint FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);

	// Read the Vertex Shader code from the file
	std::string VertexShaderCode;
	std::ifstream VertexShaderStream(vertex_file_path, std::ios::in);
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





bool OpenGL_Renderer::init( int argc, char** argv)
{
    
    if( !glfwInit() )
    {
            fprintf( stderr, "Failed to initialize GLFW\n" );
                return false;
    }


    //glfwWindowHint(GLFW_SAMPLES, 4);//Antialiasing
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    //glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    //glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    // Open a window and create its OpenGL context
    //GLFWwindow* window; // this variable is global for simplicity)
    window = glfwCreateWindow( 1024, 768, "Tutorial 01", NULL, NULL);
    if( window == NULL ){
        fprintf( stderr, "Failed to open GLFW window.\n" );
        glfwTerminate();
        return false;
    }
    glfwMakeContextCurrent(window); // Initialize GLEW
    glewExperimental=true; // Needed in core profile
    if (glewInit() != GLEW_OK) {
        fprintf(stderr, "Failed to initialize GLEW\n");
        return false;
    }

    

    programID = LoadShaders("./vertex_shader.vs","./fragment_shader.fs");
    glUseProgram(programID);


	MatrixID = glGetUniformLocation(programID, "MVP");
	Projection = glm::perspective(glm::radians(45.0f), 4.0f / 3.0f, 0.1f, 100.0f);
	View       = glm::lookAt(
                        glm::vec3(4,3,3), // Camera is at (4,3,3), in World Space
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

    int nparticles = 3;
    static const GLfloat g_particle_x[] = { 
		0.2f,
		0.6f,
		0.0f,
    };
    static const GLfloat g_particle_y[] = { 
		0.2f,
		0.0f,
		0.1f,
    };
    static const GLfloat g_particle_z[] = { 
		0.2f,
		0.5f,
		0.5f,
    };
    glPointSize(8.);
	
    
    //GLuint box_VAO;
	glGenVertexArrays(1, &box_VAO);
	glBindVertexArray(box_VAO);

    
	//GLuint box_points_VBO;
	glGenBuffers(1, &box_points_VBO);
	glBindBuffer(GL_ARRAY_BUFFER, box_points_VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data), g_vertex_buffer_data, GL_STATIC_DRAW);


	glGenVertexArrays(1, &particles_VAO);
	glBindVertexArray(particles_VAO);

	glGenBuffers(1, &particles_VBO);
	glBindBuffer(GL_ARRAY_BUFFER, particles_VBO);
	glBufferData(GL_ARRAY_BUFFER, nparticles*3*sizeof(GLfloat), (void*) NULL, GL_DYNAMIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER,
                    0,
                    nparticles*sizeof(GLfloat),
                    g_particle_x);
	glBufferSubData(GL_ARRAY_BUFFER,
                    nparticles*sizeof(GLfloat),
                    nparticles*sizeof(GLfloat),
                    g_particle_y);
	glBufferSubData(GL_ARRAY_BUFFER,
                    2*nparticles*sizeof(GLfloat),
                    nparticles*sizeof(GLfloat),
                    g_particle_z);
	//glBufferData(GL_ARRAY_BUFFER, 0, nparticles*sizeof(GLfloat), GL_STATIC_DRAW);

}



//We don't wanna modify
void OpenGL_Renderer::draw()
{

    glClear( GL_COLOR_BUFFER_BIT );

    //draw_particles(render_info.n_particles, render_info.x, render_info.y, render_info.z);
    printf("draw\n");
    draw_box();

	// Swap buffers

}

void OpenGL_Renderer::draw_box()
{

    
		// Clear the screen
		glClear( GL_COLOR_BUFFER_BIT );

		// Use our shader
		glUseProgram(programID);

		// Send our transformation to the currently bound shader, 
		// in the "MVP" uniform
		glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);
        

        glEnableVertexAttribArray(0);
		glBindVertexArray(box_VAO);
		glBindBuffer(GL_ARRAY_BUFFER, box_points_VBO);
		glVertexAttribPointer(
			0,                  // attribute. No particular reason for 0, but must match the layout in the shader.
			3,                  // size
			GL_FLOAT,           // type
			GL_FALSE,           // normalized?
			0,                  // stride
			(void*)0            // array buffer offset
		);

		// Draw the triangle !
		glDrawArrays(GL_LINES, 0, 2*12); // 3 indices starting at 0 -> 1 triangle
        glDisableVertexAttribArray(0);

        glEnableVertexAttribArray(0);
        //Draw the points
		glBindVertexArray(particles_VAO);
		glBindBuffer(GL_ARRAY_BUFFER, particles_VBO);

        int nparticles = 3;
		glVertexAttribPointer(
			0,                  // attribute. No particular reason for 0, but must match the layout in the shader.
			3,                  // size
			GL_FLOAT,           // type
			GL_FALSE,           // normalized?
			sizeof(GLfloat)*nparticles,// stride
			(void*)0            // array buffer offset
		);

		// Draw the triangle !
		glDrawArrays(GL_POINTS, 0, 3); // 3 indices starting at 0 -> 1 triangle


		//glDisableVertexAttribArray(0);

   
	glfwSwapBuffers(window);
}

void OpenGL_Renderer::draw_particles(   size_t n_particles,
                                        const float *x,
                                        const float *y,
                                        const float *z
                                    )
{
    for( size_t i= 0; i < n_particles; i++){
        float posX = x[i];
        float posY = y[i];
        float posZ = z[i];
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

