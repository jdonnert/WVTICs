#include "../globals.h"

void setup_Constant_Density()
{
    Problem.Boxsize[0] = 1;
    Problem.Boxsize[1] = 1;
    Problem.Boxsize[2] = 1;

    sprintf ( Problem.Name, "IC_Constant_Density" );

    Density_Func_Ptr = &Constant_Density;
}

float Constant_Density ( const int ipart )
{
    double volume = Problem.Boxsize[0] * Problem.Boxsize[1] * Problem.Boxsize[2];
    double mass = Param.Npart * Problem.Mpart;

    return mass / volume;
}
