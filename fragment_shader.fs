#version 330 core

uniform vec3 fragment_color;

in vec3 frag_color;
in vec3 frag_normal;

out vec3 color;
void main(){



  color = frag_color;
  //color = vec3(1.,0.,0.);
}
