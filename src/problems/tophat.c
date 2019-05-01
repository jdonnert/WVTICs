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

    Density_Func_Ptr = &TopHat_Density;
}

float TopHat_Density ( const int ipart , const double bias )
{
    float x = P[ipart].Pos[0] / Problem.Boxsize[0];

    double ret; // store density return value for bias correction
    // double bias = 0.17273624122661307;  // for 1e4
    // double bias = 0.132452809143783;        // for 1e6

    const double halfstep = DENSITY_STEP * RHO_MEAN;

    const double rho_max = RHO_MEAN + halfstep;
    const double rho_min = RHO_MEAN - halfstep;

    if ( x <= 0.1 || x > 0.9 ) {
        ret = rho_min;
        ret += (ret - RHO_MEAN) * bias;
        return ret;
    } else if ( x > 0.4 && x <= 0.6 ) {
        ret = rho_max;
        ret += (ret - RHO_MEAN) * bias;
        return ret;
    } else if ( x > 0.6 ) {
        ret = rho_max - ( rho_max - rho_min ) * ( x - 0.6 ) / 0.3;
        ret += (ret - RHO_MEAN) * bias;
        return ret;
    }

    ret = rho_min + ( rho_max - rho_min ) * ( x - 0.1 ) / 0.3;
    ret += (ret - RHO_MEAN) * bias;
    return ret;
}
