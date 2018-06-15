#include "../globals.h"

// Adjust all these to your liking
void setup_User_Density ( const int subflag )
{
    Problem.Boxsize[0] = 300.0;
    Problem.Boxsize[1] = 200.0;
    Problem.Boxsize[2] = 1.0;

    Problem.Periodic[0] = false;

    sprintf ( Problem.Name, "IC_User" );

    const double rho0 = 0.0027747626;
    const double xmin = 3.0;
    const double af = 0.78655113;
    const double xne0 = 0.13;
    Problem.Rho_Max = rho0 * xne0 * exp ( -af * log ( xmin ) );
    printf ( "RhoMax = %g\n", Problem.Rho_Max );

    Density_Func_Ptr = &User_Density;
    U_Func_Ptr = &User_U;
    Velocity_Func_Ptr = &User_Velocity;
    Magnetic_Field_Func_Ptr = &User_Bfld;
    PostProcessing_Func_Ptr = &User_PostProcessing;
}

float User_Density ( const int ipart )
{
    const double rho0 = 0.0027747626;
    const double xmin = 3.0;
    const double af = 0.78655113;
    const double xne0 = 0.13;

    const double x = P[ipart].Pos[0];

    return rho0 * xne0 * exp ( -af * log ( x + xmin ) );
}

float User_U ( const int ipart )
{
    return 1.0f;
}

void User_Velocity ( const int ipart, float out[3] )
{
    out[0] = out[1] = out[2] = 0;
}

void User_Bfld ( const int ipart, float out[3] )
{
    out[0] = out[1] = out[2] = 0;
}

void User_PostProcessing()
{
}
