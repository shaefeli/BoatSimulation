#include <iostream>
#include <vector>
#include <fstream>
using namespace std;

class Maya_Interface{
public:
    string fileStart;
    string fileEnd;
    string ids;
    Maya_Interface(size_t nr_particles);
    void writeToMaya(size_t frameNr, float* x ,float* y, float* z, size_t nr_particles);

};