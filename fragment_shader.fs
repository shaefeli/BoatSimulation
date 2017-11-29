#version 330 core

uniform vec3 fragment_color;

in vec3 frag_color;

out vec3 color;
void main(){
  //color = fragment_color;
  color = frag_color;
}
