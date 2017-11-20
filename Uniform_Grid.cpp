#include "Uniform_Grid.h"
#include <cstdio>
#include <cmath>
#include <vector>

Uniform_Grid::Uniform_Grid(float min_x, float min_y, float min_z,
             float max_x, float max_y, float max_z,
             float cell_x,float cell_y,float cell_z):
             min_x(min_x),   min_y(min_y),   min_z(min_z),
             max_x(max_x),   max_y(max_y),   max_z(max_z),
             cell_x(cell_x), cell_y(cell_y), cell_z(cell_z)
{
    n_cells_x = ceil((max_x-min_x)/cell_x);
    n_cells_y = ceil((max_y-min_y)/cell_y);
    n_cells_z = ceil((max_z-min_z)/cell_z);
    cells = std::vector<std::vector<size_t>>(n_cells_x*n_cells_y*n_cells_z,
                                                std::vector<size_t>(0));
}

size_t Uniform_Grid::unroll_grid_position( size_t i, size_t j, size_t k )
{
    //printf("cell sizes: %lu %lu %lu\n", n_cells_x, n_cells_y, n_cells_z);
    return n_cells_y*n_cells_z*i + n_cells_z*j + k;
}

void Uniform_Grid::build(   const float *xs,
                            const float *ys,
                            const float *zs,
                            size_t n_particles )
{

    for( size_t i = 0; i < n_particles; i++ ) {
        float x,y,z;
        x = xs[i];
        y = ys[i];
        z = zs[i];

        //Normalize
        x = (x-min_x)/(max_x-min_x);
        y = (y-min_y)/(max_y-min_y);
        z = (z-min_z)/(max_z-min_z);

        //Move from [0,1] -> "grid space"
        //Clamp it to the correct range, so if it's outside we consider it to be on the edge
        //box. Should work fine, but better to clamp positions on the other code
        size_t gi,gj,gk;
        gi = std::max(size_t(0),std::min(n_cells_x-1,size_t(floor(x*n_cells_x))));
        gj = std::max(size_t(0),std::min(n_cells_y-1,size_t(floor(y*n_cells_y))));
        gk = std::max(size_t(0),std::min(n_cells_z-1,size_t(floor(z*n_cells_z))));

        size_t position = unroll_grid_position(gi,gj,gk);
        //printf("position: %lu\n", position);
        //printf("length: %lu\n", cells.size());
        cells[position].push_back(i);
    }
}

void Uniform_Grid::query_neighbors(float x, float y, float z,
                        std::vector<size_t> &cells)
{

    //Normalize
    x = (x-min_x)/(max_x-min_x);
    y = (y-min_y)/(max_y-min_y);
    z = (z-min_z)/(max_z-min_z);

    //Move from [0,1] -> "grid space"
    size_t gi,gj,gk;
    gi = floor(x*n_cells_x);
    gj = floor(y*n_cells_y);
    gk = floor(z*n_cells_z);
    
    cells = std::vector<size_t>(0);
    for( int i = -1; i < 2; i++ ) {
        for( int j = -1; j < 2; j++ ) {
            for( int k = -1; k < 2; k++ ) {
                size_t ti,tj,tk;
                ti = gi-i;
                tj = gj-j;
                tk = gk-k;
                if( (ti >= 0 and ti < n_cells_x) and
                    (tj >= 0 and tj < n_cells_y) and
                    (tk >= 0 and tk < n_cells_z) ) {

                    cells.push_back( unroll_grid_position(ti,tj,tk) );
                }
            }
        }
    }
    
    
}

void Uniform_Grid::query_cell(float x, float y, float z, size_t &i, size_t &j, size_t &k)
{
        x = std::max(min_x,x);
        y = std::max(min_y,y);
        z = std::max(min_z,z);
        //Normalize
        x = (x-min_x)/(max_x-min_x);
        y = (y-min_y)/(max_y-min_y);
        z = (z-min_z)/(max_z-min_z);

        //Move from [0,1] -> "grid space"
        i = floor(x*n_cells_x);
        j = floor(y*n_cells_y);
        k = floor(z*n_cells_z);
}



