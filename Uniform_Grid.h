#include <vector>
#include <cstdlib>
#include <cstddef>

#ifndef _PBS_UNIFORM_GRID_
#define _PBS_UNIFORM_GRID_

class Uniform_Grid
{
    public:
        //Best way do this non-manually, can look into improvements
        //std::vector<std::vector<size_t>> cells;
        size_t *cell_info;  //Data contained by the cells
        size_t **cells;     //Pointer to the cells with the indices
        size_t *cell_size;  //Number of particles in a cell
        size_t n_cells;     //Number of cells
        size_t n_particles;
        
        bool filled;//If true there is particle information in the cells

        size_t n_cells_x, n_cells_y, n_cells_z;
        float min_x, min_y, min_z;
        float max_x, max_y, max_z;
        float cell_x, cell_y, cell_z;
        size_t unroll_grid_position( size_t i, size_t j, size_t k );
        Uniform_Grid(float min_x, float min_y, float min_z,
                     float max_x, float max_y, float max_z,
                     float cell_x,float cell_y,float cell_z);
        ~Uniform_Grid();

        //Builds all the corresponding structures
        void build(const float *x, const float *y, const float *z, size_t n_particles);
        
        //Returns a vector with the position of the neinghboring cells in the vector
        void query_neighbors(float x, float y, float z, std::vector<size_t> &cells);
        
        //Returns the indices of the particle x,y,z in the grid, as i,j,k
        void query_cell(float x, float y, float z, size_t &i, size_t &j, size_t &k);
        
        //cleans all the information contained
        //Used internally mostly
        void clean_up();
};

#endif
