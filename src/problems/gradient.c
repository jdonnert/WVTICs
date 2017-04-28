#include "../globals.h"

#define DENSITY_STEP 0.5
#define RHO_MEAN 1.0

void setup_Gradient_Density()
{
    sprintf ( Problem.Name, "IC_GradientDensity" );

    Density_Func_Ptr = &Gradient_Density;
}

float Gradient_Density ( const int ipart )
{
    float x = P[ipart].Pos[0] / Problem.Boxsize[0];

    const double halfstep = DENSITY_STEP * RHO_MEAN;

    const double rho_max = RHO_MEAN + halfstep;
    const double rho_min = RHO_MEAN - halfstep;

    if ( x <= 0.25 ) {
        return rho_min;
    }
    if ( x >= 0.75 ) {
        return rho_max;
    }
    return rho_min + ( rho_max - rho_min ) * ( x - 0.25 ) / 0.5;
}
