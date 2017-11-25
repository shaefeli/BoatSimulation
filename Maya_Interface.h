#include <iostream>
using namespace std;

class Maya_Interface{
public:
    stringstream fileStart;
    stringstream fileEnd;
    static void init(size_t nr_particles);
    static void writeToMaya(size_t frameNr, float* interleaved_positions, size_t nr_particles);

};