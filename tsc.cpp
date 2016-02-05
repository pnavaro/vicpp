#include "tsc.hpp"

TSC::TSC(const int npar)
{
    px = (double *) malloc(sizeof(double)*npar);
    py = (double *) malloc(sizeof(double)*npar);
    pz = (double *) malloc(sizeof(double)*npar);

    vx = (double *) malloc(sizeof(double)*npar);
    vy = (double *) malloc(sizeof(double)*npar);
    vz = (double *) malloc(sizeof(double)*npar);

    wx = (double *) malloc(sizeof(double)*npar);
    wy = (double *) malloc(sizeof(double)*npar);
    wz = (double *) malloc(sizeof(double)*npar);

    sx = (double *) malloc(sizeof(double)*npar);
    sy = (double *) malloc(sizeof(double)*npar);
    sz = (double *) malloc(sizeof(double)*npar);

    dpx = (double *) malloc(sizeof(double)*npar);
    dpy = (double *) malloc(sizeof(double)*npar);
    dpz = (double *) malloc(sizeof(double)*npar);

    vol = (double *) malloc(sizeof(double)*npar);

    ipx = (int *) malloc(sizeof(int)*npar);
    ipy = (int *) malloc(sizeof(int)*npar);
    ipz = (int *) malloc(sizeof(int)*npar);

    // calcul des coefficients du schema tsc
    // tsc (i) = a(i)*x^2 + b(i)*x  + c(i)   : i = 0, 1, 2


    double ctsc1 = -0.625;
    double ctsc2 = 0.75;
    for(int i=1;i<=3;i++)
    {
        double alpha = i-2.;
        aco[i] = ctsc1*alpha*alpha+ctsc2;
        bco[i] = 0.5*alpha;
        cco[i] = 1.5*alpha*alpha-1.;
    }

}


void TSC::proj(const Grid &G)
//Calcul de W sur le maillage
{

    memset(G.wx,0.,G.nx*G.ny*G.nz*sizeof(double));
    memset(G.wy,0.,G.nx*G.ny*G.nz*sizeof(double));
    memset(G.wz,0.,G.nx*G.ny*G.nz*sizeof(double));

    for(int kk=0; kk<3; kk++)
    {
        for(int i=0;i<npar;i++)
            lt[i] = ipz[i] + kk-1;

        for(int jj=0; jj<3; jj++)
        {
            for(int i=0;i<npar;i++)
                kt[i] = ipy[i] + jj-1;

            for(int ii=0; ii<3; ii++)
            {
                for(int i=0;i<npar;i++)
                    jt[i] = ipx[i] + ii-1;

                for(int i=0;i<npar;i++)
                {

                    double sgr =   1.0
                            * (aco[ii]+dpx[i]*(bco[ii]+cco[ii]*dpx[i]))
                            * (aco[jj]+dpy[i]*(bco[jj]+cco[jj]*dpy[i]))
                            * (aco[kk]+dpz[i]*(bco[kk]+cco[kk]*dpz[i]));

                    int l = lt[i]+G.nz*(kt[i]+G.ny*jt[i]);

                    G.wx[l] += wx[i]*sgr;
                    G.wy[l] += wy[i]*sgr;
                    G.wz[l] += wz[i]*sgr;

                    G.vol[l] += vol[i]*sgr;

                }
            }
        }
    }

    double dxyz = G.dx*G.dy*G.dz;

    for (int l=0;l<G.nx*G.ny*G.nz;l++)
    {
        G.wx[l] /= dxyz;
        G.wy[l] /= dxyz;
        G.wz[l] /= dxyz;
    }

}

void TSC::interp(const Grid & G)
//Calcul de V et S sur les particules
{

    memset(sx,0.,npar*sizeof(double));
    memset(sy,0.,npar*sizeof(double));
    memset(sz,0.,npar*sizeof(double));


    for(int kk=0; kk<3; kk++)
    {
        for(int i=0;i<npar;i++) lt[i] = ipz[i] + kk-1;
        for(int jj=0; jj<3; jj++)
        {
            for(int i=0;i<npar;i++) kt[i] = ipy[i] + jj-1;
            for(int ii=0; ii<3; ii++)
            {
                for(int i=0;i<npar;i++) jt[i] = ipx[i] + ii-1;

                for(int i=0;i<npar;i++)
                {

                    double sgr = (aco[ii]+dpx[i]*(bco[ii]+cco[ii]*dpx[i]))
                            * (aco[jj]+dpy[i]*(bco[jj]+cco[jj]*dpy[i]))
                            * (aco[kk]+dpz[i]*(bco[kk]+cco[kk]*dpz[i]));

                    int l = lt[i]+G.nz*(kt[i]+G.ny*jt[i]);

                    vx[i] += G.vx[l]*sgr;
                    vy[i] += G.vy[l]*sgr;
                    vz[i] += G.vz[l]*sgr;

                    sx[i] += G.sx[l]*sgr;
                    sy[i] += G.sy[l]*sgr;
                    sz[i] += G.sz[l]*sgr;

                }
            }
        }
    }

    double dxyz = G.dx*G.dy*G.dz;

    for (int l=0;l<G.nx*G.ny*G.nz;l++)
    {
        sx[l] *=  dxyz;
        sy[l] *=  dxyz;
        sz[l] *=  dxyz;
    }

}

void TSC::spot(const Grid & G)
//Reperage des particules ip (numero du noeud)
//dp : position de la particule dans la maille
{
    for(int i=0;i<npar;++i)
    {
        ipx[i] = int((px[i]-G.xmin)/G.dx);
        ipy[i] = int((py[i]-G.ymin)/G.dy);
        ipz[i] = int((pz[i]-G.zmin)/G.dz);

        dpx[i] = (px[i]-ipx[i]*G.dx-G.xmin)/G.dx;
        dpy[i] = (py[i]-ipy[i]*G.dy-G.ymin)/G.dy;
        dpz[i] = (pz[i]-ipz[i]*G.dz-G.zmin)/G.dz;

        if( ipx[i] < 1 && ipx[i] > G.nx-2 )
        {
            cout << "erreur "<< endl;
            cout << " i     " << i;
            cout << " ipx   " << ipx[i];
            cout << " x     " << px[i];
            cout << " dx    " << G.dx;
            cout << " xmin  " << G.xmin << endl;
            exit(1);
        }
        if( ipy[i] < 1 && ipy[i] > G.ny-2 )
        {
            cout << "erreur "<< endl;
            cout << " i    " << i;
            cout << " ipy  " << ipy[i];
            cout << " y    " << py[i];
            cout << " dy   " << G.dy;
            cout << " ymin " << G.ymin << endl;
        }
        if( ipz[i] < 1 && ipz[i] > G.nz-2 )
        {

            cout << "erreur "<< endl;
            cout << " i    " << i;
            cout << " ipz  " << ipz[i];
            cout << " z    " << pz[i];
            cout << " dz   " << G.dz;
            cout << " zmin " << G.zmin << endl;
        }
    }

}
