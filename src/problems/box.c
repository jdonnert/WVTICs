#include "../globals.h"

void setup_Box()
{
    Problem.Boxsize[0] = 1.0;
    Problem.Boxsize[1] = 1.0;
    Problem.Boxsize[2] = 0.1;

    sprintf ( Problem.Name, "Box" );

    const double rho1 = 4.0;
    const double rho2 = 1.0;
    const double rho_mean = 0.5 * ( rho1 + rho2 );

    Problem.Rho_Max = rho1;
    Problem.Mpart = rho_mean * ( Problem.Boxsize[0] * Problem.Boxsize[1] * Problem.Boxsize[2] ) / Param.Npart;

    Density_Func_Ptr = &Box_Density;
    U_Func_Ptr = &Box_U;
    Velocity_Func_Ptr = &Box_Velocity;
}

bool isInnerBox ( const int ipart )
{
    const double x = P[ipart].Pos[0];
    const double y = P[ipart].Pos[1];
    const double z = P[ipart].Pos[2];

    if ( x <= Problem.Boxsize[0] / 2 && y <= Problem.Boxsize[0] / 2 && Problem.Boxsize[2] / 2 ) {
        return true;
    } else {
        return false;
    }
}

float Box_Density ( const int ipart )
{
    if ( isInnerBox ( ipart ) ) {
        return 4.0;
    } else {
        return 1.0;
    }
}


/* The next step is setting up the velocity profile for the Box, following for example Hopkins 2015 or Hu 2014 */


void Box_Velocity ( const int ipart, float out[3] )
{
    out[0] = 142.3;
    out[1] = -31.3;
    out[2] = 0.0;
}

/* We set up the internal energy via the pressure profile and ideal equation of state */

float Box_U ( const int ipart )
{

    const double gamma = 5.0 / 3.0;
    const double rho1 = 4.0;
    const double rho2 = 1.0;
    const double pressure = 2.5;

    if ( isInnerBox ( ipart ) ) {
        return pressure / ( gamma - 1 ) * rho1;
    } else {
        return pressure / ( gamma - 1 ) * rho2;
    }
}

/* Just a note at the end, the Box test is very tricky for any Code, because of several reasons. In a grid Code you suffer from advection errors (in principle bulk motion along the grid). For SPH th test is difficult, because of the jumping density
at the edge of the inner Box and the bad behaviour of SPH at contact discontinuities. Further we use a glass which could be problematic in terms of regularity. We would advise to use this Test with caution and probablt setting the ICs up by using a regular grid
structure and run the simulation with a moving mesh or the advanced MFM and MFV methods presented in Gaburov 2011 */

