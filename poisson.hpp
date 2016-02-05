#ifndef POISSON_HPP
#define POISSON_HPP
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <iostream>
#include "grid.hpp"

class Poisson
{
private:
    /* In place computation     */
    double * c;
    double * d;
    double  cx, cy, cz;
    fftw_plan plan_forward, plan_backward;
    double * u; /* psi */
    double * f; /* rhs */

public:  

    Poisson(const Grid & G);
    void setU( double * v, int size); // Set field
    void setF( double * v, int size); // Set RHS
    double * getU();                  // Get Solution
    void solve(const Grid & G);
    double error( const Grid & G, const double * v);
    void finalize();

    ~Poisson(); //Destructor
};

#endif // POISSON_HPP

