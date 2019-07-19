#include "../globals.h"

void setup_Evrard_Collapse()
{
    Problem.Boxsize[0] = 10.0;
    Problem.Boxsize[1] = 10.0;
    Problem.Boxsize[2] = 10.0;

    sprintf ( Problem.Name, "IC_Evrard_Collapse" );

    Problem.Rho_Max = 10; // This value is empiric

    Density_Func_Ptr = &Evrard_Collapse_Density;
    U_Func_Ptr = &Evrard_Collapse_U;
    Velocity_Func_Ptr = &Evrard_Collapse_Velocity;

}


/* At first we set up a constant density in the Box */
float Evrard_Collapse_Density ( const int ipart , const double bias )
{
    double const x = P[ipart].Pos[0] - Problem.Boxsize[0] * 0.5;
    double const y = P[ipart].Pos[1] - Problem.Boxsize[1] * 0.5;
    double const z = P[ipart].Pos[2] - Problem.Boxsize[2] * 0.5;

    double Radius = sqrt ( x * x + y * y + z * z );
    double epsilon = 0.01;

    if ( Radius < 1 ) {
        return 1.0 / ( 2 * pi * ( Radius + epsilon ) );
    } else {
        return 0.001;
    }

}

void Evrard_Collapse_Velocity ( const int ipart, float out[3] )
{
    out[0] = 0.0;
    out[1] = 0.0;
    out[2] = 0.0;
}

/* We set up the internal energy via the pressure profile and ideal equation of state */

float Evrard_Collapse_U ( const int ipart )
{
    return 0.05;
}
