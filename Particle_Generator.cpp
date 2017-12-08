#include "Particle_Generator.h"
#include <iostream>
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"







void load_model_data()
{

    std::string inputfile = "Models/watercraftPack_001.obj";
    std::string mtlfolder = "Models/";
    
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
      
    std::string err;
    bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, inputfile.c_str(), mtlfolder.c_str());
      
    if (!err.empty()) { // `err` may contain warning message.
      std::cerr << err << std::endl;
    }

    // Loop over shapes
    for (size_t s = 0; s < shapes.size(); s++) {
      // Loop over faces(polygon)
      size_t index_offset = 0;
      for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
        int fv = shapes[s].mesh.num_face_vertices[f];

        // Loop over vertices in the face.
        for (size_t v = 0; v < fv; v++) {
          // access to vertex
          tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
          tinyobj::real_t vx = attrib.vertices[3*idx.vertex_index+0];
          tinyobj::real_t vy = attrib.vertices[3*idx.vertex_index+1];
          tinyobj::real_t vz = attrib.vertices[3*idx.vertex_index+2];
          tinyobj::real_t nx = attrib.normals[3*idx.normal_index+0];
          tinyobj::real_t ny = attrib.normals[3*idx.normal_index+1];
          tinyobj::real_t nz = attrib.normals[3*idx.normal_index+2];
          tinyobj::real_t tx = attrib.texcoords[2*idx.texcoord_index+0];
          tinyobj::real_t ty = attrib.texcoords[2*idx.texcoord_index+1];
          // Optional: vertex colors
          // tinyobj::real_t red = attrib.colors[3*idx.vertex_index+0];
          // tinyobj::real_t green = attrib.colors[3*idx.vertex_index+1];
          // tinyobj::real_t blue = attrib.colors[3*idx.vertex_index+2];
        }
        index_offset += fv;

        // per-face material
        shapes[s].mesh.material_ids[f];
      }
    }

}






//Object centered around the origin
void generate_particle_cube(
                            float length,  //Length of cube side
                            float h,        //Distance between particles
                            std::vector<float> &xv,
                            std::vector<float> &yv,
                            std::vector<float> &zv,
                            size_t &n_particles 
                            )
{
    float xmin = -length;
    float ymin = -length;
    float zmin = -length;
    float xmax =  length;
    float ymax =  length;
    float zmax =  length;
    
    size_t x_particles = (xmax-xmin)/h;
    size_t y_particles = (ymax-ymin)/h;
    size_t z_particles = (zmax-zmin)/h;
    xv = std::vector<float>(x_particles*y_particles*z_particles);
    yv = std::vector<float>(x_particles*y_particles*z_particles);
    zv = std::vector<float>(x_particles*y_particles*z_particles);
    for( size_t i = 0; i < x_particles; i++ ) {
        for( size_t j = 0; j < y_particles; j++ ) {
            for( size_t k = 0; k < z_particles; k++ ) {
                xv[i*y_particles*z_particles + j*z_particles + k] = i*h+xmin;
                yv[i*y_particles*z_particles + j*z_particles + k] = j*h+ymin;
                zv[i*y_particles*z_particles + j*z_particles + k] = k*h+zmin;
            }
        }
    }
    n_particles = x_particles*y_particles*z_particles;
}

