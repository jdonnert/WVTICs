#include "../globals.h"

#define DENSITY_STEP 0.5
#define RHO_MEAN 1.0

void setup_TopHat_Density()
{
    Problem.Boxsize[0] = 1;
    Problem.Boxsize[1] = 0.5;
    Problem.Boxsize[2] = 0.1;

    sprintf ( Problem.Name, "IC_TopHat" );

    Problem.Rho_Max = 1.5;
    Problem.Mpart = RHO_MEAN * ( Problem.Boxsize[0] * Problem.Boxsize[1] * Problem.Boxsize[2] ) / Param.Npart;

    Density_Func_Ptr = &TopHat_Density;
}

float TopHat_Density ( const int ipart )
{
    float x = P[ipart].Pos[0] / Problem.Boxsize[0];

    const double halfstep = DENSITY_STEP * RHO_MEAN;

    const double rho_max = RHO_MEAN + halfstep;
    const double rho_min = RHO_MEAN - halfstep;

    if ( x <= 0.1 || x > 0.9 ) {
        return rho_min;
    } else if ( x > 0.4 && x <= 0.6 ) {
        return rho_max;
    } else if ( x > 0.6 ) {
        return rho_max - ( rho_max - rho_min ) * ( x - 0.6 ) / 0.3;
    }

    return rho_min + ( rho_max - rho_min ) * ( x - 0.1 ) / 0.3;
}
