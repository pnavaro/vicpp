#include "poisson.hpp"

Poisson::Poisson(const Grid & G)
{
    u = (double *) malloc(G.nx*G.ny*G.nz*sizeof(double));
    f = (double *) malloc(G.nx*G.ny*G.nz*sizeof(double));

    int i, j, k;
    double vx, vy, vz;

    int MX = G.nx-2;
    int MY = G.ny-2;
    int MZ = G.nz-2;

    c = (double *) malloc(MX*MY*MZ*sizeof(double));
    d = (double *) malloc(MX*MY*MZ*sizeof(double));

    cx = 1. / (G.dx*G.dx);
    cy = 1. / (G.dy*G.dy);
    cz = 1. / (G.dz*G.dz);

    for(i = 0; i < MX; i++)
        for(j = 0; j < MY; j++)
            for(k = 0; k < MZ; k++)
            {
                vx = 1.-cos((i+1)*M_PI/(MX+1.));
                vy = 1.-cos((j+1)*M_PI/(MY+1.));
                vz = 1.-cos((k+1)*M_PI/(MZ+1.));
                d[k+MZ*(j+MY*i)] = 2.*(cx*vx+cy*vy+cz*vz);
            }

//    if ( fftw_init_threads() == 0 )
//    {
//        printf("\n problem with threads \n");
//    }
//    else
//    {
//        fftw_plan_with_nthreads(2);
//    }

    plan_forward  = fftw_plan_r2r_3d(MX,MY,MZ,c,c,
                                     FFTW_RODFT00,FFTW_RODFT00,FFTW_RODFT00,
                                     FFTW_ESTIMATE);

    plan_backward = fftw_plan_r2r_3d(MX,MY,MZ,c,c,
                                     FFTW_RODFT00,FFTW_RODFT00,FFTW_RODFT00,
                                     FFTW_ESTIMATE);

    printf(" \n Initialisation du solveur de Poisson \n ");
}

void Poisson::solve(const Grid & G)
{
    int i, j, k;

    int NY = G.ny;
    int NZ = G.nz;

    int MX = G.nx-2;
    int MY = G.ny-2;
    int MZ = G.nz-2;

    /* Boundary conditions */
    for(i=0; i < MX; i++)
        for(j=0; j < MY; j++)
            for(k=0; k < MZ; k++)
            {

                int l = k+MZ*(j+MY*i);

                c[l] = f[(k+1)+NZ*((j+1)+NY*(i+1))];

                if (i==0 )   c[l] += cx * u[(k+1)+NZ*((j+1)+NY*(0))];
                if (i==MX-1) c[l] += cx * u[(k+1)+NZ*((j+1)+NY*(NY-1))];
                if (j==0 )   c[l] += cy * u[(k+1)+NZ*((0)+NY*(i+1))];
                if (j==MY-1) c[l] += cy * u[(k+1)+NZ*((NY-1)+NY*(i+1))];
                if (k==0 )   c[l] += cz * u[(0)+NZ*((j+1)+NY*(i+1))];
                if (k==MZ-1) c[l] += cz * u[(NZ-1)+NZ*((j+1)+NY*(i+1))];
            }

    fftw_execute( plan_forward );

    for(i=0; i < MX; i++)
        for(j=0; j < MY; j++)
            for(k=0; k < MZ; k++)
            {
                c[k+MZ*(j+MY*i)] = c[k+MZ*(j+MY*i)] / d[k+MZ*(j+MY*i)];
            }

    fftw_execute( plan_backward );

    // normalisation
    for(i=0; i < MX; i++)
        for(j=0; j < MY; j++)
            for(k=0; k < MZ; k++)
            {
                u[(k+1)+NZ*((j+1)+NY*(i+1))] = c[k+MZ*(j+MY*i)]/((2*MX+2)*(2*MY+2)*(2*MZ+2));
            }
}

void Poisson::setU( double * v, int size)
{
    for (int l=0; l<size; l++)
        this->u[l] = v[l];
}

void Poisson::setF( double * v, int s)
{
    for (int l=0; l<s; l++)
    {
        f[l] = v[l];
    }
}

double * Poisson::getU()
{
    return this->u;
}

double Poisson::error ( const Grid & G, const double * v)
{
double t;
double e;
for(int i = 0; i < G.nx; i++)
    for(int j = 0; j < G.ny; j++)
       for(int k = 0; k < G.nz; k++)
       {
           int l = k+G.nz*(j+G.ny*i);
           t = fabs(this->u[l]-v[l]);
           if ( t > e ) e = t;
       }
 return e;
}


Poisson::~Poisson()
{
    finalize();
}

void Poisson::finalize()
{
    //free memory
    fftw_destroy_plan( plan_forward );
    fftw_destroy_plan( plan_backward );
//    fftw_cleanup_threads();
    free(c);
    free(d);
}

