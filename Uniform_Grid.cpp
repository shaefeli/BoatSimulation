#include "Uniform_Grid.h"
#include <cassert>
#include <queue>
#include <cstdio>
#include <iostream>
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
    //cells = std::vector<std::vector<size_t>>(n_cells_x*n_cells_y*n_cells_z,
                                                //std::vector<size_t>(0));
    cells = (size_t **)malloc(sizeof(size_t *)*n_cells_x*n_cells_y*n_cells_z);
    cell_size = (size_t *)malloc(sizeof(size_t)*n_cells_x*n_cells_y*n_cells_z);
    n_cells = n_cells_x*n_cells_y*n_cells_z;

    std::cerr<<"grid size:"<<n_cells_x<<"x"<<n_cells_y<<"x"<<n_cells_z<<std::endl;

    filled = false;
}
Uniform_Grid::~Uniform_Grid()
{
    //std::cerr<<"destroying grid"<<std::endl;
    if(filled) clean_up();
    free(cells);
    free(cell_size);
    //std::cerr<<"grid destroyed"<<std::endl;
}

size_t Uniform_Grid::unroll_grid_position( size_t i, size_t j, size_t k )
{
    //printf("cell sizes: %lu %lu %lu\n", n_cells_x, n_cells_y, n_cells_z);
    assert( 0 <= i and i < n_cells_x );
    assert( 0 <= j and j < n_cells_y );
    assert( 0 <= k and k < n_cells_z );
    return n_cells_y*n_cells_z*i + n_cells_z*j + k;
}

void Uniform_Grid::build(   const float *xs,
                            const float *ys,
                            const float *zs,
                            size_t n_particles )
{
    //Free the memory from the previous run
    if(filled) clean_up();

    std::vector<std::queue<size_t>> cell_queues(n_cells);
    //std::cerr<<"Pushing particles"<<std::endl;
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

        gi = floor(x*n_cells_x);
        gj = floor(y*n_cells_y);
        gk = floor(z*n_cells_z);

        //Clamp it to a correct value
        gi = std::max(size_t(0),std::min(n_cells_x-1,gi));
        gj = std::max(size_t(0),std::min(n_cells_y-1,gj));
        gk = std::max(size_t(0),std::min(n_cells_z-1,gk));

        size_t position = unroll_grid_position(gi,gj,gk);

        cell_queues[position].push(i);
    }
    //std::cerr<<"Finished pushing particles"<<std::endl;

    for( size_t i = 0; i < n_cells; i++ ) {

        cell_size[i] = cell_queues[i].size();
        cells[i] = (size_t *)malloc(sizeof(size_t)*cell_size[i]);

        size_t *cell = cells[i];

        size_t j = 0;
        while( not cell_queues[i].empty() ) {
            //cells[i][j] = cell_queues[i].front();
            cell[j] = cell_queues[i].front();

            //Check that values are correct
            assert( cell[j] < n_particles );
            //std::cerr<<"part:"<<cell_queues[i].front()<<std::endl;
            cell_queues[i].pop();
            j++;
        }
        assert( j == cell_size[i] );
        //std::cerr<<"All elements added: "<<j<<" "<<cell_size[i]<<std::endl;
        //std::cerr<<"Finished adding particles into the array"<<std::endl;
    }
    filled = true;
    //std::cerr<<"Filled positions"<<std::endl;
}

void Uniform_Grid::clean_up()
{
    //std::cerr<<"clean up"<<std::endl;
    for( int i = 0; i < n_cells; i++ ) {
        free(cells[i]);
        cell_size[i] = 0;
    }
    filled = false;
    //std::cerr<<"clean"<<std::endl;
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

    //clamp
    gi = std::max(size_t(0),std::min(size_t(n_cells_x-1),gi));
    gj = std::max(size_t(0),std::min(size_t(n_cells_y-1),gj));
    gk = std::max(size_t(0),std::min(size_t(n_cells_z-1),gk));
    
    cells = std::vector<size_t>(0);
    for( int i = -1; i < 2; i++ ) {
        for( int j = -1; j < 2; j++ ) {
            for( int k = -1; k < 2; k++ ) {
                size_t ti,tj,tk;
                ti = gi + i;
                tj = gj + j;
                tk = gk + k;
                if( (ti >= 0 and ti < (n_cells_x)) and
                    (tj >= 0 and tj < (n_cells_y)) and
                    (tk >= 0 and tk < (n_cells_z)) ) {
                    //std::cerr<<"neighbors:"<<ti<<" "<<tj<<" "<<tk<<std::endl;
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



