#include "vortexring.hpp"

VortexRing::VortexRing(int npar)
    : TSC(npar)
{
    //  Tore de grand rayon rtore et de petit rayon rdisc.
    //  l'axe de revolution du tore est (oz). le plan meridien
    //  est : x = 0
    //
    //  le tore est la surface balayee par le disque de centre
    //  (0,rmoy,0) et de rayon ray, lorsqu'il effectue des
    //	rotations d'angle phi autour de l'axe de revolution.
    //

    printf("\n circulation par anneau       : gama  = %lf12.3",gama);
    printf("\n rayon du tore                : rtore = %lf12.3",rtore);
    printf("\n rayon de la section          : rdisc = %lf12.3",rdisc);
    printf("\n nb de rayons de la section   : nray  = %d12   ",nray);

    double *rf = nullptr;
    double *zf = nullptr;
    double *of = nullptr;
    double *sf = nullptr;
    double *h  = nullptr; //distance moyenne entre deux particules

    //distribution de la vorticite initiale dans la section
    nfil = this->distri_vort_init (rf, zf, sf, of, h );

    //dlmax = 000.0; dlmin = 100.0

    //npano = int(2.*pi*rmoy/h)
    //npar = nfil * npano
    //write(*,'("nombre de particules par filaments : npano =",i12  )') npano
    //write(*,'("nombre de particules total         : npar  =",i12  )') npar
    //allocate(px(npar),py(npar),pz(npar),gamp(npar),volp(npar))
    //allocate(wx(npar),wy(npar),wz(npar))

    //m = 0
    //do k = 1, nfil     !loop over filaments

    //   dlong = 2. * pi * rf(k) / float(npano)
    //   dlmax = amax1(dlmax, dlong)
    //   dlmin = amin1(dlmin, dlong)

    //   dphi = 2. * pi / float(npano)
    //   do j = 1, npano
    //      phi     = -pi + (j-1) * dphi
    //      m       = m+1
    //      px(m)   = rf(k) * cos(phi)
    //      py(m)   = rf(k) * sin(phi)
    //      pz(m)   = zf(k)
    //      gamp(m) = of(k)
    //      volp(m) = dlong * sf(k)
    //      !wx(m)   = -sin(phi)*of(k)*dlong
    //      !wy(m)   =  cos(phi)*of(k)*dlong
    //      !wz(m)   = 0.0
    //   end do

    //end do

    //!*** calcul du poids initial de chaque particule

    //j2=0
    //do i = 1, nfil
    //   j1 = j2+1
    //   j2 = j1+npano-1
    //   do j = j1, j2
    //      if (j == j1) then
    //         i1=j2  ; i2=j1+1
    //      else if (j == j2) then
    //         i1=j2-1; i2=j1
    //      else
    //         i1=j-1; i2=j+1
    //      end if
    //      wx(j) = gamp(j)*(.5*(px(i2)-px(i1)))
    //      wy(j) = gamp(j)*(.5*(py(i2)-py(i1)))
    //      wz(j) = 0.0
    //   end do
    //end do

    //! *** volume de l'anneau ***

    //volum = 2.0 * pi * rmoy * pi * ray * ray

    //write(*,*) 'volumes theorique - numerique', volum, sum(volp)
    //write(*,*) 'circulation totale', sum(gamp)
    //ut = ( log(8.*rmoy/ray) - 0.25 ) / ( 4.*pi*rmoy )
    //write(*,"('vitesse theorique de translation du tore : ',f7.5,/)")ut

    //!*** pas d'espace : distance interparticulaire

    //hmax = amax1( h, dlmax)
    //hmin = amin1( h, dlmin)

    //write(*,"('distance entre 2 particules dans un filament')")
    //write(*,"('distance maxi ',f12.6, ' mini ',f12.6 )") dlmax, dlmin
    //write(*,"('distance entre 2 particules dans le tore')")
    //write(*,"('distance maxi ',f12.6, ' mini ',f12.6 )") hmax, hmin
    //write(*,*)

    //deallocate(rf,zf,of,sf)

}




//end subroutine init_vortexring

int VortexRing::distri_vort_init (double * rf, double * zf, double * sf,
                                  double * of, double * h )
{
    // dr    : pas d'espace dans la direction radiale
    // nray  : nb de pts ds la direction radiale.
    // dray  : rayon de la particule placee au centre
    // rdisc : rayon de la section
    // gama  : circulation totale
    // surf  : surface de la section
    // nsec  : nombre de points dans la section

    double r1, s1, r2, s2, r, dss, teta, q;
    int  nsec0;
    bool gauss = false;

    double dr      = rdisc / ( nray + 0.5 );
    double dray    = 0.5 * dr;	//rayon de la section centrale
    double surf    = M_PI * rdisc * rdisc;
    nsec    = 6;
    double dteta   = 2.0 * M_PI / nsec;
    double dazimax = dr * dteta;


    //estimate the particle number
    int k = 1;
    for(int i=0;i<nray;i++)
        k += i*nsec;

    rf = (double *) malloc(k*sizeof(double));
    sf = (double *) malloc(k*sizeof(double));
    of = (double *) malloc(k*sizeof(double));

    rf[0] = 0.0;
    zf[0] = 0.0;
    sf[0] = M_PI * dray * dray;
    if (gauss)
        of[0] = gama * ( 1.-exp(-(dray/rdisc)*(dray/rdisc))) ;
    else
        of[0] = gama * sf[0] / surf	;

    //write(57,1000) rf(1), zf(1), of(1)

    r1   = dray;
    s1   = M_PI * r1*r1;

    nsec0 = nsec;
    nsec  = 0;

    for(int i=0;i<nray;i++)
    {

        nsec    = nsec + nsec0;
        dteta   = 2.0 * M_PI / nsec;
        r       = ( i+1 ) * dr;
        if (r*dteta>dazimax)
            dazimax= r * dteta ;

        r2      = r + 0.5 * dr;
        s2      = M_PI * r2*r2;
        dss     = s2 - s1;
        s1      = s2;

        for(int j=0;j<nsec; j++)
        {

            k ++;

            teta    = (j+1) * dteta;
            rf[k] = r * cos( teta );
            zf[k] = r * sin( teta );
            sf[k] = dss /  nsec ;	//surface d'un secteur

            if (gauss)
            {
                q = 0.5 * (exp(-(r1/rdisc)*(r1/rdisc))-exp(-(r2/rdisc)*(r2/rdisc)));
                of[k]= gama * dteta / M_PI * q;
            }
            else
            {
                of[k] = gama * sf[k] / surf;
            }

            printf("\n %lf %lf %lf ",rf[k], zf[k], of[k]);

        }

        r1 = r2;
        int kd = k - nsec + 1;
        printf("\n %lf %lf %lf ", rf[kd], zf[kd], of[kd]);

    }

    return k;

    *h = sqrt( dr*dr + dazimax * dazimax );

    //printf('nb de pts sur la section : %d)', nf);
    //printf('----------------------------------------------');
    //printf('distance entre 2 secteurs voisins:            ');
    //printf('.....dans la direction radiale   : %lf ',    dr);
    //printf('.....dans la direction azimutale : %lf ', dazimax);
    //printf('------------------------------');
    // printf('surface theorique - pratique : %lf %lf ', surf,sum(sf));

    //if (gauss) then
    //   write(*,*)'gama theorique - pratique :',gama*(1.-exp(-1.)),' - ',sum(of)
    //else
    //   write(*,*)'gama theorique - pratique :',gama, ' - ',sum(of)
    //end if

}
