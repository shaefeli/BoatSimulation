#include "Uniform_Grid.h"
#include <algorithm>
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
             //max_x(max_x+cell_x),   max_y(max_y),   max_z(max_z),
             cell_x(cell_x), cell_y(cell_y), cell_z(cell_z)
{
    n_particles = 0;
    cell_info = NULL;
    n_cells_x = ceil((max_x-min_x)/cell_x);
    n_cells_y = ceil((max_y-min_y)/cell_y);
    n_cells_z = ceil((max_z-min_z)/cell_z);
    //Allow exceeding bounds through the positive side
    //Just in case the correspondence is wrong
    this->max_x = n_cells_x*cell_x+min_x;
    this->max_y = n_cells_y*cell_y+min_y;
    this->max_z = n_cells_z*cell_z+min_z;
    
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
    free(cell_info);
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
    //Move into clean up maybe
    if( this->n_particles != n_particles ) {
        this->n_particles = n_particles;
        if( cell_info != NULL ) free(cell_info);
        cell_info = (size_t *)malloc(sizeof(size_t)*n_particles);
    }

    //Free the memory from the previous run
    if(filled) clean_up();
    
    //First position: cell index
    //Second position: particle index
    std::vector<std::pair<size_t,size_t>> particle_index(n_particles);
    size_t gi,gj,gk;
    for( size_t i = 0; i < n_particles; i++ ) {
        query_cell(xs[i],ys[i],zs[i],gi,gj,gk);
        particle_index[i].first  = unroll_grid_position(gi,gj,gk);
        particle_index[i].second = i;
    }
    std::sort(particle_index.begin(),particle_index.end(),
                            []( const std::pair<size_t,size_t> &v1,
                                const std::pair<size_t,size_t> &v2) -> bool {
                                return v1.first < v2.first;
                            });

    //copy sorted info
    for( int i = 0; i < n_particles; i++ ) {
        cell_info[i] = particle_index[i].second;
    }

    std::vector<std::pair<size_t,size_t>>::iterator low,high;
    low = particle_index.begin();
    for( int i = 0; i < n_cells; i++ ) {
        //std::cerr<<"pos:"<<i<<std::endl;
        high = upper_bound(particle_index.begin(),particle_index.end(),
                            std::pair<size_t,size_t>(i,0),
                            []( const std::pair<size_t,size_t> &v1,
                                const std::pair<size_t,size_t> &v2) -> bool {
                                return v1.first < v2.first;
                            });

        size_t n_elements = high - low;
        //cells[i] = (size_t *)malloc(sizeof(size_t)*n_elements);
        cells[i] = &(cell_info[low-particle_index.begin()]);
        cell_size[i] = n_elements;

        //std::cerr<<"n_elements:"<<n_elements<<std::endl;
        /*
        for( int j = 0; j < n_elements; j++ ) {
            //std::cerr<<"vals "<<i<<":"<<low->first<<" "<<low->second<<std::endl;
            assert( low->first == i );
            cells[low->first][j] = low->second;
            low++;
        }
        */
        low = high;
    }
    //std::cerr<<"done"<<std::endl;
    filled = true;
    //std::cerr<<"Filled positions"<<std::endl;
}

void Uniform_Grid::clean_up()
{
    /*
    //std::cerr<<"clean up"<<std::endl;
    for( int i = 0; i < n_cells; i++ ) {
        free(cells[i]);
        cell_size[i] = 0;
    }
    filled = false;
    //std::cerr<<"clean"<<std::endl;
    */
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
    //If they get outside the bounds, they get assigned to the closest in-bounds
    gi = size_t(std::max(0.,floor(x*n_cells_x)));
    gj = size_t(std::max(0.,floor(y*n_cells_y)));
    gk = size_t(std::max(0.,floor(z*n_cells_z)));

    //clamp
    gi = std::min(size_t(n_cells_x-1),gi);
    gj = std::min(size_t(n_cells_y-1),gj);
    gk = std::min(size_t(n_cells_z-1),gk);
    
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
                    std::cerr<<"neighbors:"<<ti<<" "<<tj<<" "<<tk<<std::endl;
                    cells.push_back( unroll_grid_position(ti,tj,tk) );
                }
            }
        }
    }
    std::cerr<<std::endl;
    
    
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
        i = size_t(std::max(0.,floor(x*n_cells_x)));
        j = size_t(std::max(0.,floor(y*n_cells_y)));
        k = size_t(std::max(0.,floor(z*n_cells_z)));

        i = std::min(size_t(n_cells_x-1),i);
        j = std::min(size_t(n_cells_y-1),j);
        k = std::min(size_t(n_cells_z-1),k);
}



