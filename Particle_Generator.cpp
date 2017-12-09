#include "Particle_Generator.h"
#include <iostream>
#include <limits>
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"


#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>

#include <thinks/poissonDiskSampling.hpp>

//typedef CGAL::Simple_cartesian<double>               Kernel;
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef Kernel::Point_3                              Point_3;
typedef Kernel::Vector_3                             Vector_3;
typedef CGAL::Polyhedron_3<Kernel>                   Polyhedron_3;
typedef CGAL::Nef_polyhedron_3<Kernel>               Nef_polyhedron_3;

typedef CGAL::Surface_mesh<Point_3> Surface_mesh;

typedef Nef_polyhedron_3::Vertex_const_handle Vertex_const_handle;
typedef Nef_polyhedron_3::Halfedge_const_handle Halfedge_const_handle;
typedef Nef_polyhedron_3::Halffacet_const_handle Halffacet_const_handle;
typedef Nef_polyhedron_3::Volume_const_handle Volume_const_handle;
typedef Nef_polyhedron_3::Volume_const_iterator Volume_const_iterator;
typedef Nef_polyhedron_3::Object_handle Object_handle;



template <typename T, std::size_t N>
class Vec
{
public:
    typedef T value_type;
    static const std::size_t size = N;
    Vec() {}
    T& operator[](std::size_t i) { return _data[i]; }
    const T& operator[](std::size_t i) const { return _data[i]; }
private:
    T _data[N];
};



void load_model_data(   float h,
                        std::vector<float> &xv,
                        std::vector<float> &yv,
                        std::vector<float> &zv,
                        size_t &n_particles 
        )
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

    std::vector<Polyhedron_3> Ps(shapes.size());

    
    float minx = std::numeric_limits<float>::max();
    float miny = std::numeric_limits<float>::max();
    float minz = std::numeric_limits<float>::max();

    float maxx = std::numeric_limits<float>::min();
    float maxy = std::numeric_limits<float>::min();
    float maxz = std::numeric_limits<float>::min();

    float avgx = 0.;
    float avgy = 0.;
    float avgz = 0.;

    for( size_t i = 0; i < attrib.vertices.size()/3; i++ ) {
        tinyobj::real_t vx = attrib.vertices[3*i+0];
        tinyobj::real_t vy = attrib.vertices[3*i+1];
        tinyobj::real_t vz = attrib.vertices[3*i+2];

        avgx += vx;
        avgy += vy;
        avgz += vz;

        minx = std::min(minx,vx); 
        miny = std::min(miny,vy);
        minz = std::min(minz,vz);

        maxx = std::max(maxx,vx);
        maxy = std::max(maxy,vy);
        maxz = std::max(maxz,vz);
    }
    avgx /= attrib.vertices.size()/3;
    avgy /= attrib.vertices.size()/3;
    avgz /= attrib.vertices.size()/3;


    //Center to 0,0
    maxx -= avgx;
    maxy -= avgy;
    maxz -= avgz;

    minx -= avgx;
    miny -= avgy;
    minz -= avgz;

   float scale = 0.25;


    // Loop over shapes
    for (size_t s = 0; s < shapes.size(); s++) { //Only work with last one for now :I
      std::vector<Point_3> points;
      // Loop over faces(polygon)
      size_t index_offset = 0;
      for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
        int fv = shapes[s].mesh.num_face_vertices[f];

        // Loop over vertices in the face.
        for (size_t v = 0; v < fv; v++) {
          // access to vertex
          tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
          tinyobj::real_t vx = attrib.vertices[3*idx.vertex_index+0] - avgx;
          tinyobj::real_t vy = attrib.vertices[3*idx.vertex_index+1] - avgy;
          tinyobj::real_t vz = attrib.vertices[3*idx.vertex_index+2] - avgz;
          vx *= scale;
          vy *= scale;
          vz *= scale;


          //Only want information about the positions of the vertices to create the polygons...
          //tinyobj::real_t nx = attrib.normals[3*idx.normal_index+0];
          //tinyobj::real_t ny = attrib.normals[3*idx.normal_index+1];
          //tinyobj::real_t nz = attrib.normals[3*idx.normal_index+2];
          //tinyobj::real_t tx = attrib.texcoords[2*idx.texcoord_index+0];
          //tinyobj::real_t ty = attrib.texcoords[2*idx.texcoord_index+1];
          points.push_back(Point_3(vx,vy,vz));

          // Optional: vertex colors
          // tinyobj::real_t red = attrib.colors[3*idx.vertex_index+0];
          // tinyobj::real_t green = attrib.colors[3*idx.vertex_index+1];
          // tinyobj::real_t blue = attrib.colors[3*idx.vertex_index+2];
        }
        index_offset += fv;

        // per-face material
        shapes[s].mesh.material_ids[f];
      }
      CGAL::convex_hull_3(points.begin(), points.end(), Ps[s]);
      std::cout << "The convex hull contains " << Ps[s].size_of_vertices() << " vertices" << std::endl;
    }
    maxx *= scale;
    maxy *= scale;
    maxz *= scale;
    minx *= scale;
    miny *= scale;
    minz *= scale;
    avgx *= scale;
    avgy *= scale;
    avgz *= scale;


    //std::vector<Nef_polyhedron_3> NPs(Ps.size());
    Nef_polyhedron_3 NP;
    for( int i = 0; i < Ps.size(); i++ ) NP += Nef_polyhedron_3(Ps[i]);

    std::cout<<avgx<<std::endl;
    std::cout<<avgy<<std::endl;
    std::cout<<avgz<<std::endl;
    std::cout<<std::endl;
    std::cout<<minx<<std::endl;
    std::cout<<miny<<std::endl;
    std::cout<<minz<<std::endl;
    std::cout<<std::endl;
    std::cout<<maxx<<std::endl;
    std::cout<<maxy<<std::endl;
    std::cout<<maxz<<std::endl;



    //Poisson sample the bounding box



    /*
    Vertex_const_handle v;
    Halfedge_const_handle e;
    Halffacet_const_handle f;
    Volume_const_handle c;

    Object_handle o = NP.locate(Point_3(0.5,0.5,0.5));
    if(CGAL::assign(v,o)) {
      std::cout << "Locating vertex" << std::endl;
    }
    else if(CGAL::assign(e,o)) {
      std::cout << "Locating edge" << std::endl;
    }
    else if(CGAL::assign(f,o)) {
      std::cout << "Locating facet" << std::endl;
    }
    else if(CGAL::assign(c,o)) {
      std::cout << "Locating volume" << std::endl;
      //for( auto vi = NP.volumes_begin(); vi != NP.volumes_end(); ++vi ) {
        //std::cout<<(vi->mark())<<std::endl;
      //}
      std::cout<<c->mark()<<std::endl;
    }
    */
    Vec<float,3> x_min, x_max;
    x_min[0] = minx;
    x_min[1] = miny;
    x_min[2] = minz;

    x_max[0] = maxx;
    x_max[1] = maxy;
    x_max[2] = maxz;
    uint32_t max_sample_attempts = 30;
    uint32_t seed = 1981;
    std::vector<Vec<float,3>> bounding_box_samples 
                    = thinks::poissonDiskSampling(  h,
                                                    x_min,
                                                    x_max,
                                                    max_sample_attempts,
                                                    seed);

   
   std::vector<float> interior_points_x;
   std::vector<float> interior_points_y;
   std::vector<float> interior_points_z;
   size_t n_interior_points = 0;
   for( size_t i = 0; i < bounding_box_samples.size(); i++ ) {
        Point_3 p(  bounding_box_samples[i][0],
                    bounding_box_samples[i][1],
                    bounding_box_samples[i][2]);
        Vertex_const_handle v;
        Halfedge_const_handle e;
        Halffacet_const_handle f;
        Volume_const_handle c;

        Object_handle o = NP.locate(p);
        if(CGAL::assign(v,o)) {
          std::cout << "Locating vertex" << std::endl;
        }
        else if(CGAL::assign(e,o)) {
          std::cout << "Locating edge" << std::endl;
        }
        else if(CGAL::assign(f,o)) {
          std::cout << "Locating facet" << std::endl;
        }
        else if(CGAL::assign(c,o)) {
          //std::cout << "Locating volume" << std::endl;
          std::cout<<c->mark();
          if( c->mark() ) {
            interior_points_x.push_back(bounding_box_samples[i][0]);
            interior_points_y.push_back(bounding_box_samples[i][1]);
            interior_points_z.push_back(bounding_box_samples[i][2]);
            n_interior_points++;
          }
        }
        

   }
   xv = interior_points_x;
   yv = interior_points_y;
   zv = interior_points_z;
   n_particles = n_interior_points;
   std::cout<<std::endl;



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

