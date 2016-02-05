#ifndef XDMF_HPP
#define XDMF_HPP
#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>
#include "grid.hpp"
#include "poisson.hpp"

#define H5FILE_NAME "field.h5"

class Xdmf
{
public:
    Xdmf(const Grid &G);
    void write_hdf5_mesh(const Grid &G);
    void write_hdf5_data(const Grid &G, double * u, double * f);
};

#endif // XDMF_HPP
