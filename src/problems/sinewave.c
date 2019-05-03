#include "../globals.h"

#define ORDER 1
#define RHO_MEAN 1.0

void setup_SineWave_Density()
{
    sprintf ( Problem.Name, "IC_SineWave" );

    Problem.Boxsize[0] = 1;
    Problem.Boxsize[1] = 0.75;
    Problem.Boxsize[2] = 0.75;

    Problem.Rho_Max = 1.5;

    Density_Func_Ptr = &SineWave_Density;
}

float SineWave_Density ( const int ipart , const double bias )
{
    float x = P[ipart].Pos[0] / Problem.Boxsize[0];

    // float bias = 0.17440653993991542;    // for 1e4
    // float bias = 0.14955494792461688;    // for 1e6
    double ret = RHO_MEAN * ( 1 + 0.5 * sin ( 2 * pi * ORDER * x ) );
    ret += (ret - RHO_MEAN) * bias;

    return ret;

}
