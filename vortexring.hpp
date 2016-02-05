#ifndef VORTEXRING_HPP
#define VORTEXRING_HPP
#include <math.h>
#include "tsc.hpp"

class VortexRing : public TSC
{
public:
    VortexRing(int npar);

private:
    int nfil  ;//  nombre de filaments
    int npano ;//  nombre segments par filament
    int nray  ;//  nombre petits rayons
    int nsec  ;//  nombre de secteurs du disque

    double gama   ;//  circulation totale
    double rtore  ;//  rayon de la section
    double rdisc  ;//  grand rayon
    double * gamp ;
    double ut;

    int distri_vort_init (double * rf, double * zf, double * sf,
                          double * of, double * h );

};

#endif // VORTEXRING_HPP
