#include "../globals.h"

void setup_Gresho_Vortex()
{
    Problem.Boxsize[0] = 1.0;
    Problem.Boxsize[1] = 1.0;
    Problem.Boxsize[2] = 0.1;

    sprintf ( Problem.Name, "IC_Gresho" );

    const double rho = 1.0;

    Problem.Rho_Max = rho;

    Density_Func_Ptr = &Gresho_Vortex_Density;
    U_Func_Ptr = &Gresho_Vortex_U;
    Velocity_Func_Ptr = &Gresho_Vortex_Velocity;
}

/* At first we set up a constant density in the Box */
float Gresho_Vortex_Density ( const int ipart )
{
    return 1.0;
}

/* Now we are writting a function to calculate the angle in polar coordinates*/

float Gresho_Vortex_Phi ( const int ipart )
{

    if ( P[ipart].Pos[0] > 0 ) {
        return atan ( P[ipart].Pos[1] / P[ipart].Pos[0] );
    } else if ( P[ipart].Pos[0] == 0 ) {
        return P[ipart].Pos[1] / abs ( P[ipart].Pos[1] ) * 0.5 * pi;
    } else if ( P[ipart].Pos[0] < 0 && P[ipart].Pos[1] >= 0 ) {
        return atan ( P[ipart].Pos[1] / P[ipart].Pos[0] ) + pi;
    } else {
        return atan ( P[ipart].Pos[1] / P[ipart].Pos[0] ) - pi;
    }
}

/* The next step is setting up the velocity profile for the Vortex, following for example Hopkins 2015 or Hu 2014 */
/* It is just the easy case with Machnumber of 1, the more difficult is done later */

void Gresho_Vortex_Velocity ( const int ipart, float out[3] )
{
    double x = P[ipart].Pos[0] - Problem.Boxsize[0] * 0.5;
    double y = P[ipart].Pos[1] - Problem.Boxsize[1] * 0.5;
    double Radius = sqrt ( x * x + y * y );

    if ( Radius < 0.2 ) {
        out[0] = -5 * Radius * sin ( Gresho_Vortex_Phi ( ipart ) );
        out[1] =  5 * Radius * cos ( Gresho_Vortex_Phi ( ipart ) );
        out[2] = 0.;
    } else if ( Radius < 0.4 ) {
        out[0] = - ( 2 - 5 * Radius ) * sin ( Gresho_Vortex_Phi ( ipart ) );
        out[1] =  ( 2 - 5 * Radius ) * cos ( Gresho_Vortex_Phi ( ipart ) );
        out[2] = 0.;
    } else {
        out[0] = 0.;
        out[1] = 0.;
        out[2] = 0.;
    }
}

/* We set up the internal energy via the pressure profile and ideal equation of state */

float Gresho_Vortex_U ( const int ipart )
{
    const double gamma = 5.0 / 3.0;
    const double rho = 1.0;
    double x = P[ipart].Pos[0] - Problem.Boxsize[0] * 0.5;
    double y = P[ipart].Pos[1] - Problem.Boxsize[1] * 0.5;
    double Radius = sqrt ( x * x + y * y );

    if ( Radius < 0.2 ) {
        return ( 5 + 12 * Radius * Radius ) / ( gamma - 1 ) / rho;
    } else if (  Radius < 0.4 ) {
        return ( 9 + 12 * Radius * Radius ) / ( gamma - 1 ) / rho;
    } else {
        return ( 3 + 4 * log ( 2 ) ) / ( gamma - 1 ) / rho;
    }
}

/* Just a note at the end, the gresho vortex is in general a twodimensional testcase. We set it up in three dimesnions, because in a common simulation we
will only use three dimensions. So you can check the third dimension, but if your code works proper it should remain zero everywhere */

