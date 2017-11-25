#include "OpenGL_Renderer.h"
#include "Maya_Interface.h"
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
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

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
	color_uniform = glGetUniformLocation(programID, "fragment_color");

	Projection = glm::perspective(glm::radians(45.0f), 4.0f / 3.0f, 0.1f, 100.0f);
	View       = glm::lookAt(
                        glm::vec3(2,2,2), // Camera is at (4,3,3), in World Space
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

    glPointSize(8.);
	
    
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
    //free(interleaved_data);//optimize this
    
    glBindVertexArray(0);

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

}

void OpenGL_Renderer::draw_box()
{
        glUniform3f(color_uniform, 1., 0., 0.);
		// Clear the screen
		glBindVertexArray(box_VAO);
        glEnableVertexArrayAttrib(box_VAO,0);
		glBindBuffer(GL_ARRAY_BUFFER, box_points_VBO);
		glVertexAttribPointer(
			0,                  // attribute
			3,                  // size
			GL_FLOAT,           // type
			GL_FALSE,           // normalized?
			0,                  // stride
			(void*)0            // array buffer offset
		);

		// Draw the triangle !
		glDrawArrays(GL_LINES, 0, 2*12); // 3 indices starting at 0 -> 1 triangle
        glDisableVertexAttribArray(0);
		glBindVertexArray(0);
}

void OpenGL_Renderer::draw_particles( )
{
        glUniform3f(color_uniform, 0., 0., 1.);

        //Draw the points
        glEnableVertexArrayAttrib(particles_VAO,0);
		glBindVertexArray(particles_VAO);
		glBindBuffer(GL_ARRAY_BUFFER, particles_VBO);

        size_t nparticles = render_info.n_particles;

        //float *interleaved_data = (float *)malloc(3*nparticles*sizeof(float));
        for( int i = 0; i < nparticles; i++ ) {
            interleaved_buffer[i*3 + 0] = render_info.x[i];
            interleaved_buffer[i*3 + 1] = render_info.y[i];
            interleaved_buffer[i*3 + 2] = render_info.z[i];
        }
     size_t frame = 0;
    Maya_Interface::writeToMaya(frame,interleaved_buffer,nparticles);
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

    glDrawArrays(GL_POINTS, 0, nparticles*3);

    glDisableVertexAttribArray(0);


}


void OpenGL_Renderer::reshape(int x, int y)
{
    if (y == 0 || x == 0) return;  //Nothing is visible then, so return
    //Set a new projection matrix
    //Angle of view:40 degrees
    //Near clipping plane distance: 0.5
    //Far clipping plane distance: 20.0
     
}

