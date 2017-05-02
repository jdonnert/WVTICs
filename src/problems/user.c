#include "../globals.h"

// Adjust all these to your liking
void setup_User_Density ( const int subflag )
{
    Problem.Boxsize[0] = 1.0;
    Problem.Boxsize[1] = 1.0;
    Problem.Boxsize[2] = 1.0;

    sprintf ( Problem.Name, "IC_User" );

    Problem.Rho_Max = 1.0;
    double rho_mean = 1.0;
    Problem.Mpart = rho_mean * ( Problem.Boxsize[0] * Problem.Boxsize[1] * Problem.Boxsize[2] ) / Param.Npart;

    Density_Func_Ptr = &User_Density;
}

float User_Density ( const int ipart )
{
    //float x = P[ipart].Pos[0] / Problem.Boxsize[0];
    //float y = P[ipart].Pos[1] / Problem.Boxsize[0];
    //float z = P[ipart].Pos[2] / Problem.Boxsize[0];

    float rho = 1.0f;

    return rho;
}