#version 330 core

// Input vertex data, different for all executions of this shader.
layout(location = 0) in vec3 vertexPosition_modelspace;
layout(location = 1) in vec3 point_color;
layout(location = 2) in vec3 normal;

out vec3 frag_color;
out vec3 frag_normal;

// Values that stay constant for the whole mesh.
uniform mat4 MVP;

uniform vec3 fragment_color;

void main(){
    vec3 color = vec3(1.,0.,0.);

    frag_normal = normal;
    //frag_color = point_color;
	// Output position of the vertex, in clip space : MVP * position
	gl_Position =  MVP * vec4(vertexPosition_modelspace,1);

    color = fragment_color;

    float an = clamp(dot(normal,normalize(vec3(1-gl_Position.x,0.2-gl_Position.y,0.-gl_Position.z))),0.,1.);
    frag_color = 0.2*color+0.8*color*an;

}

