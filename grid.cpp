#include "grid.hpp"
#define I(i,j,k) (k)+nz*((j)+ny*(i))

Grid::Grid()
{
    Grid G(0.,1.,0.,1.,0.,1.,34,34,34);
}

Grid::Grid( double xmin, double xmax,
            double ymin, double ymax,
            double zmin, double zmax,
            int nx, int ny, int nz)

{
    printf("\n nx, ny, nz = %d %d %d  ", nx, ny, nz);
    printf("\n xmin, xmax = %lf %lf   ", xmin, xmax);
    printf("\n ymin, ymax = %lf %lf   ", ymin, ymax);
    printf("\n zmin, zmax = %lf %lf \n", zmin, zmax);

    this->xmin = xmin;
    this->xmax = xmax;
    this->ymin = ymin;
    this->ymax = ymax;
    this->zmin = zmin;
    this->zmax = zmax;

    this->nx = nx;
    this->ny = ny;
    this->nz = nz;

    dx = (xmax-xmin)/(nx-1.);
    dy = (ymax-ymin)/(ny-1.);
    dz = (zmax-zmin)/(nz-1.);

    x = (double *) malloc(nx*ny*nz*sizeof(double));
    y = (double *) malloc(nx*ny*nz*sizeof(double));
    z = (double *) malloc(nx*ny*nz*sizeof(double));

    for(int i = 0; i < nx; i++)
        for(int j = 0; j < ny; j++)
            for(int k = 0; k < nz; k++)
            {
                int l = k+nz*(j+ny*i);
                x[l] = xmin + i*dx;
                y[l] = ymin + j*dy;
                z[l] = zmin + k*dz;
            }

    vx = (double *) malloc(nx*ny*nz*sizeof(double));
    vy = (double *) malloc(nx*ny*nz*sizeof(double));
    vz = (double *) malloc(nx*ny*nz*sizeof(double));

    sx = (double *) malloc(nx*ny*nz*sizeof(double));
    sy = (double *) malloc(nx*ny*nz*sizeof(double));
    sz = (double *) malloc(nx*ny*nz*sizeof(double));

    wx = (double *) malloc(nx*ny*nz*sizeof(double));
    wy = (double *) malloc(nx*ny*nz*sizeof(double));
    wz = (double *) malloc(nx*ny*nz*sizeof(double));

}

void Grid::vitdef(double * psx, double * psy, double * psz) const
{

    double dpzdy, dpydz, dpxdz, dpzdx, dpydx, dpxdy;

    double s2dx = 0.5/dx;
    double s2dy = 0.5/dy;
    double s2dz = 0.5/dz;

    static double du[4][4];


    for(int i=1;i<nx-1;i++)
        for(int j=1;j<ny-1;j++)
            for(int k=1;k<nz-1;k++)
            {
                dpzdy = psz[I(i,j+1,k)] - psz[I(i,j-1,k)];
                dpydz = psy[I(i,j,k+1)] - psy[I(i,j,k-1)];

                dpxdz = psx[I(i,j,k+1)] - psx[I(i,j,k-1)];
                dpzdx = psz[I(i+1,j,k)] - psz[I(i-1,j,k)];

                dpydx = psy[I(i+1,j,k)] - psy[I(i-1,j,k)];
                dpxdy = psx[I(i,j+1,k)] - psx[I(i,j-1,k)];

                vx[I(i,j,k)] = s2dy * dpzdy - s2dz * dpydz;
                vy[I(i,j,k)] = s2dz * dpxdz - s2dx * dpzdx;
                vz[I(i,j,k)] = s2dx * dpydx - s2dy * dpxdy;
            }

    for(int i=1;i<nx-1;i++)
        for(int j=1;j<ny-1;j++)
            for(int k=1;k<nz-1;k++)
            {
                du[1][1] = s2dx * (vx[I(i+1,j,k)] - vx[I(i-1,j,k)]);
                du[1][2] = s2dy * (vx[I(i,j+1,k)] - vx[I(i,j-1,k)]);
                du[1][3] = s2dz * (vx[I(i,j,k+1)] - vx[I(i,j,k-1)]);

                du[2][1] = s2dx * (vy[I(i+1,j,k)] - vy[I(i-1,j,k)]);
                du[2][2] = s2dy * (vy[I(i,j+1,k)] - vy[I(i,j-1,k)]);
                du[2][3] = s2dz * (vy[I(i,j,k+1)] - vy[I(i,j,k-1)]);

                du[3][1] = s2dx * (vz[I(i+1,j,k)] - vz[I(i-1,j,k)]);
                du[3][2] = s2dy * (vz[I(i,j+1,k)] - vz[I(i,j-1,k)]);
                du[3][3] = s2dz * (vz[I(i,j,k+1)] - vz[I(i,j,k-1)]);

                // s = du^t . w

                int l = I(i,j,k);
                sx[I(i,j,k)] = du[1][1]*wx[l]+du[2][1]*wy[l]+du[3][1]*wz[l];
                sy[I(i,j,k)] = du[1][2]*wx[l]+du[2][2]*wy[l]+du[3][2]*wz[l];
                sz[I(i,j,k)] = du[1][3]*wx[l]+du[2][3]*wy[l]+du[3][3]*wz[l];

            }
}
