#include "Particle_Generator.h"
#include <iostream>
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"


#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Nef_polyhedron_3.h>

#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/convex_decomposition_3.h> 


//typedef CGAL::Simple_cartesian<double>               Kernel;
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef Kernel::Point_3                              Point_3;
typedef CGAL::Polyhedron_3<Kernel>                   Polyhedron;
typedef Polyhedron::Facet_iterator                   Facet_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_facet_circulator;

typedef CGAL::Nef_polyhedron_3<Kernel, CGAL::SNC_indexed_items> Nef_polyhedron_3;
typedef Nef_polyhedron_3::Volume_const_iterator Volume_const_iterator;



typedef Polyhedron::HalfedgeDS             HalfedgeDS;

template <class HDS>
class My_Builder : public CGAL::Modifier_base<HDS> {
public:
    tinyobj::attrib_t attrib;
    tinyobj::shape_t shape;
    std::vector<tinyobj::material_t> materials;
    My_Builder(
            tinyobj::attrib_t attrib,
            tinyobj::shape_t &shape,
            std::vector<tinyobj::material_t> &materials
            ) {
        this->attrib = attrib;
        this->shape = shape;
        this->materials = materials;
    }
    void operator()( HDS& hds) {
      // Postcondition: hds is a valid polyhedral surface.
      /*
      CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
      B.begin_surface( 3, 1, 6);
      typedef typename HDS::Vertex   Vertex;
      typedef typename Vertex::Point Point;
      B.add_vertex( Point( 0, 0, 0));
      B.add_vertex( Point( 1, 0, 0));
      B.add_vertex( Point( 0, 1, 0));
      B.begin_facet();
      B.add_vertex_to_facet( 0);
      B.add_vertex_to_facet( 1);
      B.add_vertex_to_facet( 2);
      B.end_facet();
      B.end_surface();
      */
      CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true );

      long n_vertices = shape.mesh.indices.size();
      long n_faces = shape.mesh.num_face_vertices.size();
      
      std::cout<<"n_vertices:"<<n_vertices<<std::endl;
      std::cout<<"n_faces:"<<n_faces<<std::endl;

      B.begin_surface(n_vertices,n_faces);
      //Add all vertices

      std::cout<<"shape: "<<shape.name<<std::endl;
      // Loop over faces(polygon)
      size_t index_offset = 0;
      for (size_t f = 0; f < shape.mesh.num_face_vertices.size(); f++) {

        int fv = shape.mesh.num_face_vertices[f];
        assert(fv == 3);//Only triangles
        // Loop over vertices in the face.
        for (size_t v = 0; v < fv; v++) {
          // access to vertex
          tinyobj::index_t idx = shape.mesh.indices[index_offset + v];
          tinyobj::real_t vx = attrib.vertices[3*idx.vertex_index+0];
          tinyobj::real_t vy = attrib.vertices[3*idx.vertex_index+1];
          tinyobj::real_t vz = attrib.vertices[3*idx.vertex_index+2];
          B.add_vertex(Point_3( vx, vy, vz));
        }
        B.begin_facet();
        for( int i = 0; i < 3; i++ ) {
            B.add_vertex_to_facet(index_offset++);
        }
        B.end_facet();
        //index_offset += fv;

        // per-face material
        //shapes[s].mesh.material_ids[f];
      }
      B.end_surface();

    }
};



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

    std::vector<Polyhedron> Ps(shapes.size());
    for( int i = 0; i < shapes.size(); i++ ) {
        My_Builder<HalfedgeDS> builder(attrib,shapes[i],materials);
        Ps[i].delegate(builder);
        std::cout<<"Is closed:"<<Ps[i].is_closed()<<std::endl;
    }

    //std::vector<Nef_polyhedron_3> nPs(Ps);

    /*
    CGAL::convex_decomposition_3(N);
    std::list<Polyhedron_3> convex_parts;
          
    //the first volume is the outer volume, which is 
    //ignored in the decomposition
    Volume_const_iterator ci = ++N.volumes_begin();
    for( ; ci != N.volumes_end(); ++ci) {
        if(ci->mark()) {
            Polyhedron_3 P;
            N.convert_inner_shell_to_polyhedron(ci->shells_begin(), P);
            convex_parts.push_back(P);
        }
    }
    */

    /*
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

          //Only want information about the positions of the vertices to create the polygons...
          //tinyobj::real_t nx = attrib.normals[3*idx.normal_index+0];
          //tinyobj::real_t ny = attrib.normals[3*idx.normal_index+1];
          //tinyobj::real_t nz = attrib.normals[3*idx.normal_index+2];
          //tinyobj::real_t tx = attrib.texcoords[2*idx.texcoord_index+0];
          //tinyobj::real_t ty = attrib.texcoords[2*idx.texcoord_index+1];

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
    */

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

