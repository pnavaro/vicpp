#ifndef TSC_HPP
#define TSC_HPP
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "grid.hpp"

using namespace std;

class TSC
{
public:
    TSC(const int npar);

    int npar;

    void proj( const Grid & G);
    void interp( const Grid & G);
    void spot( const Grid & G);

    void step( const double dt);

private:
    double * px;
    double * py;
    double * pz;

    double * vx;
    double * vy;
    double * vz;

    double * wx;
    double * wy;
    double * wz;

    double * sx;
    double * sy;
    double * sz;

    int * jt;
    int * kt;
    int * lt;
    int * nt;

    int * ipx;
    int * ipy;
    int * ipz;

    double * dpx;
    double * dpy;
    double * dpz;

    double aco[3];
    double bco[3];
    double cco[3];

    double * vol;


};

#endif // TSC_HPP
