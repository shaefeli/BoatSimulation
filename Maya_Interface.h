#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
using namespace std;

class Maya_Interface{
public:
    string fileStart;
    string fileEnd;
    string ids;
    int maxIter;
    int minIter;
    
    float *mob_x;
    float *mob_y;
    float *mob_z;

    float *rot_x;
    float *rot_y;
    float *rot_z;

    Maya_Interface(size_t nr_particles, int maxIt, int minIt);
    void writeToMaya(size_t frameNr, float* x ,float* y, float* z, size_t nr_particles);
    void writeBoatInfos(float x, float y, float z, float rx, float ry, float rz, int it);

    ~Maya_Interface();

};