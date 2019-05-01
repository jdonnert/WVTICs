#include "../globals.h"

#define DENSITY_STEP 0.5
#define RHO_MEAN 1.0

void setup_Sawtooth_Density()
{
    Problem.Boxsize[0] = 1;
    Problem.Boxsize[1] = 0.1;
    Problem.Boxsize[2] = 0.1;

    sprintf ( Problem.Name, "IC_Sawtooth" );

    Problem.Rho_Max = 1.5;

    Density_Func_Ptr = &Sawtooth_Density;
}

float Sawtooth_Density ( const int ipart , const double bias = 0.0)
{
    float x = P[ipart].Pos[0] / Problem.Boxsize[0];
    double ret;
    // float bias = 0.10176706576254425; // 1e4
    // float bias = 0.019724626821411465; // 1e6
    if ( x > 0.5 ) {
        x -= 0.5;
    }

    const double halfstep = DENSITY_STEP * RHO_MEAN;

    const double rho_max = RHO_MEAN + halfstep;
    const double rho_min = RHO_MEAN - halfstep;

    ret = rho_min + ( rho_max - rho_min ) * x / 0.5;
    ret += (ret - RHO_MEAN) * bias;

    return ret;
}
