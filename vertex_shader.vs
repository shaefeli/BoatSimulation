#version 330 core

// Input vertex data, different for all executions of this shader.
layout(location = 0) in vec3 vertexPosition_modelspace;
layout(location = 1) in vec3 point_color;

out vec3 frag_color;

// Values that stay constant for the whole mesh.
uniform mat4 MVP;

void main(){
    frag_color = point_color;
	// Output position of the vertex, in clip space : MVP * position
	gl_Position =  MVP * vec4(vertexPosition_modelspace,1);

}

