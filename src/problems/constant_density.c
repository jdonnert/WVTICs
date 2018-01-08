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
    return 1.0;
}
