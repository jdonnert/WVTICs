#include "../globals.h"

#define ORDER 4
#define RHO_MEAN 1.0

void setup_SineWave_Density()
{
    sprintf ( Problem.Name, "IC_SineWave" );

    Problem.Boxsize[0] = 1;
    Problem.Boxsize[1] = 0.75;
    Problem.Boxsize[2] = 0.75;

    Problem.Rho_Max = 1.5;
    Problem.Mpart = RHO_MEAN * ( Problem.Boxsize[0] * Problem.Boxsize[1] * Problem.Boxsize[2] ) / Param.Npart;

    Density_Func_Ptr = &SineWave_Density;
}

float SineWave_Density ( const int ipart )
{
    float x = P[ipart].Pos[0] / Problem.Boxsize[0];

    return RHO_MEAN * ( 1 + 0.5 * sin ( 2 * pi * ORDER * x ) );
}
