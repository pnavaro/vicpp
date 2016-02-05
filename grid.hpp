#ifndef GRID_HPP
#define GRID_HPP
#include <stdlib.h>
#include <hdf5.h>

class Grid
{
public:  //Tous les attributs devraient être privés!

    int nx;
    int ny;
    int nz;

    Grid();
    Grid(double xmin, double xmax,
         double ymin, double ymax,
         double zmin, double zmax,
         int nx, int ny, int nz);

    void vitdef(double * psx, double * psy, double * psz) const;

    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zmin;
    double zmax;

    double dx;
    double dy;
    double dz;

    double * x;    /* X coordinates  */
    double * y;    /* Y coordinates  */
    double * z;    /* Z coordinates  */

    double * vx;
    double * vy;
    double * vz;

    double * wx;
    double * wy;
    double * wz;

    double * sx;
    double * sy;
    double * sz;

    double * vol;


private:



};

#endif // GRID_HPP
