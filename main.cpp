#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#include <hdf5.h>
#include <fftw3.h>

#include "grid.hpp"
#include "poisson.hpp"
#include "xdmf.hpp"
#include "vortexring.hpp"

int
main (void)
{
    double * v;    /* Exact solution field  */

    Grid G(0.,1.,0.,1.,0.,1.,34,34,34);

    int NX = G.nx;
    int NY = G.ny;
    int NZ = G.nz;

    Poisson P(G);
    //VortexRing VR(100);
    //VR.spot(G);

    v = new double[NX*NY*NZ]; //(double *) malloc(sizeof(double)*NX*NY*NZ);

    for(int i = 0; i < NX; i++)
        for(int j = 0; j < NY; j++)
           for(int k = 0; k < NZ; k++)
           {
               int l = k+NZ*(j+NY*i);
               v[l] = cos(2*M_PI*G.x[l])
                    * sin(2*M_PI*G.y[l])
                    * cos(2*M_PI*G.z[l]);
           }

    P.setU(v,NX*NY*NZ);

    for(int i = 0; i < NX; i++)
        for(int j = 0; j < NY; j++)
           for(int k = 0; k < NZ; k++)
           {
               int l = k+NZ*(j+NY*i);
               v[l] = 12.*M_PI*M_PI * v[l];
           }


    P.setF(v,NX*NY*NZ);

    P.solve(G);

    for(int i = 0; i < NX; i++)
        for(int j = 0; j < NY; j++)
           for(int k = 0; k < NZ; k++)
           {
               int l = k+NZ*(j+NY*i);
               v[l] = cos(2*M_PI*G.x[l])
                    * sin(2*M_PI*G.y[l])
                    * cos(2*M_PI*G.z[l]);
           }

    printf( " \n error = %lf \n ", P.error(G, v) );
    //Plot results
    Xdmf XMF(G);
    XMF.write_hdf5_mesh(G);
    XMF.write_hdf5_data(G, P.getU(), v);

return 1;
}
