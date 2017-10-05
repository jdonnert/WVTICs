#include "../globals.h"

void setup_Kelvin_Helmholtz_Instability()
{
    Problem.Boxsize[0] = 256;
    Problem.Boxsize[1] = 256;
    Problem.Boxsize[2] = 16;

    sprintf ( Problem.Name, "IC_KelvinHelmholtz" );

    Problem.Rho_Max = 6.26E-8;

    Density_Func_Ptr = &Kelvin_Helmholtz_Instability_Density;
    U_Func_Ptr = &Kelvin_Helmholtz_Instability_U;
    Velocity_Func_Ptr = &Kelvin_Helmholtz_Instability_Velocity;
}

bool isOuterLayer ( const int ipart )
{
    const double y = P[ipart].Pos[1];
    const double frac = y / Problem.Boxsize[1];

    if ( frac <= 1.0 / 3.0 || frac > 2.0 / 3.0 ) {
        return true;
    } else {
        return false;
    }
}

float Kelvin_Helmholtz_Instability_Density ( const int ipart )
{
    if ( isOuterLayer ( ipart ) ) {
        return 3.13E-8;
    } else {
        return 6.26E-8;
    }
}

float Kelvin_Helmholtz_Instability_U ( const int ipart )
{
    return 101527.0;
}

void Kelvin_Helmholtz_Instability_Velocity ( const int ipart, float out[3] )
{
    if ( isOuterLayer ( ipart ) ) {
        out[0] = 40.0;
    } else {
        out[0] = -40.0;
    }

    const double deltaVy = 4.0;
    const double lambda = 128.0;
    const double x = P[ipart].Pos[0];
    const double y = P[ipart].Pos[1];

    out[1] = deltaVy * ( sin ( 2.0 * pi * ( x + lambda / 2.0 ) / lambda ) * exp ( -pow ( 10.0 * ( y - 64.0 ), 2.0 ) )
                         - sin ( 2.0 * pi * x / lambda ) * exp ( -pow ( 10.0 * ( y + 64.0 ), 2.0 ) ) );

    out[2] = 0.0;
}
