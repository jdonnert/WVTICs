#include "../globals.h"

#define DENSITY_STEP 0.5

float Gradient_Density ( const int ipart )
{
    float x = P[ipart].Pos[0] / Problem.Boxsize[0];

    double volume = Problem.Boxsize[0] * Problem.Boxsize[1] * Problem.Boxsize[2];
    double mass = Param.Npart * Problem.Mpart;

    const double rho0 = mass / volume;
    const double halfstep = DENSITY_STEP * rho0;

    const double rho_max = rho0 + halfstep;
    const double rho_min = rho0 - halfstep;

    if ( x <= 0.25 ) {
        return rho_min;
    }
    if ( x >= 0.75 ) {
        return rho_max;
    }
    return rho_min + ( rho_max - rho_min ) * ( x - 0.25 ) / 0.5;
}
