//
// Created by Sergey Ivannikov on 10/12/2017.
//

#ifndef PBS_PROJECT_SUPPORTINGSTRUCTURES_H
#define PBS_PROJECT_SUPPORTINGSTRUCTURES_H

template <typename T, std::size_t N>
class Vec
{
public:
    typedef T value_type;
    static const std::size_t size = N;
    Vec() {}
    Vec(T x, T y, T z){ _data[0] = x; _data[1] = y; _data[2] = z; }
    T& operator[](std::size_t i) { return _data[i]; }
    const T& operator[](std::size_t i) const { return _data[i]; }
private:
    T _data[N];
};


typedef struct{
    float x1,x2;
    float y1,y2;
    float z1,z2;
} BoundaryBox;

typedef struct{
    float x1,x2;
    float y1,y2;
    float z1,z2;
    float spawningRadius;
} ParticlesInitialSpawningBox;

typedef struct {
    float dt;
    float kernel_radius;
    float particle_radius;
    float rho0;
    float k;
    float mu;
    float g;
    float mass;
    float gamma; // surface tension
} SimState;

typedef struct {
    float cells_x, cells_y, cells_z;
} UniformGridSplit;


typedef struct ParticlesAllInfoStorage{                    //particles of this type have index:
    unsigned int n_liquid_particles_start;
    unsigned int n_liquid_particles;      //[n_liquid_particles_start... + n_liquid_particles]

    unsigned int n_boundary_particles_start;
    unsigned int n_boundary_particles;    //[n_boundary_particles_start... + n_boundary_particles]

    unsigned int n_mobile_particles_start;
    unsigned int n_mobile_particles;      //[n_mobile_particles_start... + n_mobile_particles]
    unsigned int n_total_particles;

    //Position
    float *x;
    float *y;
    float *z;

    // Intermediate Position
    float *x_star;
    float *y_star;
    float *z_star;

    //Velocity
    float *vx;
    float *vy;
    float *vz;

    // Intermediate Velocity
    float *vx_star;
    float *vy_star;
    float *vz_star;

    // !Total force acting
    float *Fx;
    float *Fy;
    float *Fz;

    float *Fx_p;
    float *Fy_p;
    float *Fz_p;

    // ! Curvatures
    float *nx;
    float *ny;
    float *nz;

    // densities
    float *rho;
    float *rho_star;

    // pressures
    float *p;

    float mass; // for now all particles have the same mass


    void allocateMemoryForNParticles(unsigned int N){
        x = new float[N];
        y = new float[N];
        z = new float[N];

        x_star = new float[N];
        y_star = new float[N];
        z_star = new float[N];

        vx = new float[N];
        vy = new float[N];
        vz = new float[N];

        vx_star = new float[N];
        vy_star = new float[N];
        vz_star = new float[N];

        Fx = new float[N];
        Fy = new float[N];
        Fz = new float[N];

        // pressure force
        Fx_p = new float[N];
        Fy_p = new float[N];
        Fz_p = new float[N];


        nx = new float[N];
        ny = new float[N];
        nz = new float[N];

        rho = new float[N];
        rho_star = new float[N];
        p   = new float[N];
    }

    void deallocateMemory(){
        delete [] x;        delete [] y;        delete [] z;
        delete [] x_star;   delete [] y_star;   delete [] z_star;

        delete [] vx;       delete [] vy;       delete [] vz;
        delete [] vx_star;  delete [] vy_star;  delete [] vz_star;

        delete [] Fx;   delete [] Fy;   delete [] Fz;
        delete [] Fx_p; delete [] Fy_p; delete [] Fz_p;

        delete [] nx;   delete [] ny;   delete [] nz;
        delete [] rho;
        delete [] rho_star;
        delete []p;
    }

    ~ParticlesAllInfoStorage(){
        this->deallocateMemory();
    }

} ParticlesSystemData;

#endif //PBS_PROJECT_SUPPORTINGSTRUCTURES_H

