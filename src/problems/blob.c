#include "../globals.h"

void setup_Blob()
{
    Problem.Boxsize[0] = 2000.0;
    Problem.Boxsize[1] = 2000.0;
    Problem.Boxsize[2] = 8000.0;

    sprintf ( Problem.Name, "IC_Blob" );

    const double rho = 3.13e-7; // This value is empiric

    Problem.Rho_Max = rho;
    Problem.Mpart = rho * ( Problem.Boxsize[0] * Problem.Boxsize[1] * Problem.Boxsize[2] )   / Param.Npart;

    Density_Func_Ptr = &Blob_Density;
    U_Func_Ptr = &Blob_U;
    Velocity_Func_Ptr = &Blob_Velocity;

}


/* At first we set up a constant density in the Box */
float Blob_Density ( const int ipart )
{
    double const x = P[ipart].Pos[0] - Problem.Boxsize[0] * 0.5;
    double const y = P[ipart].Pos[1] - Problem.Boxsize[1] * 0.5;
    double const z = P[ipart].Pos[2] - 3000.0;

    double Radius = sqrt ( x * x + y * y + z * z );

    if ( Radius < 197.0 ) {
        return 3.13e-7;
    } else {
        return 3.13e-8;
    }

}

void Blob_Velocity ( const int ipart, float out[3] )
{
    out[0] = 1000.0;
    out[1] = 0.0;
    out[2] = 0.0;
}

/* We set up the internal energy via the pressure profile and ideal equation of state */

float Blob_U ( const int ipart )
{
    return 0.05;
}

