#include "../globals.h"

void setup_Keplerian_Ring ()
{
    Problem.Boxsize[0] = 4.0;
    Problem.Boxsize[1] = 4.0;
    Problem.Boxsize[2] = 0.4;

    sprintf ( Problem.Name, "IC_KeplerianRing" );

    const double rho = 0.01 + ( 1 / 0.5 ) * ( 1 / 0.5 );

    Problem.Rho_Max = rho;

    Density_Func_Ptr = &Keplerian_Ring_Density;
    U_Func_Ptr = &Keplerian_Ring_U;
    Velocity_Func_Ptr = &Keplerian_Ring_Velocity;
}

/* At first we set up a constant density in the Box */
float Keplerian_Ring_Density ( const int ipart , const double bias )
{

    const double x = P[ipart].Pos[0] - Problem.Boxsize[0] * 0.5;
    const double y = P[ipart].Pos[1] - Problem.Boxsize[1] * 0.5;
    const double Radius = sqrt ( x * x + y * y );

    if ( Radius < 0.5 ) {
        return 0.01 + Radius / 0.5 / 0.5 / 0.5;
    } else if ( Radius >= 0.5 && Radius <= 2.0 ) {
        return 0.01 + 1;
    } else {
        return 0.01 + 1 / pow ( ( 1 + ( Radius - 2 ) ) / 0.1, 3 );
    }

}


/* We set up the internal energy via the pressure profile and ideal equation of state */
float Keplerian_Ring_U ( const int ipart )
{
    const double gamma = 5.0 / 3.0;
    const double Pressure = 1e-6;

    return Pressure / SphP[ipart].Rho / ( gamma - 1 );

}

void Keplerian_Ring_Velocity ( const int ipart, float out[3] )
{

    out[0] = 0.0;
    out[1] = 0.0;
    out[2] = 0.0;

}

/* Just a note at the end, the gresho vortex is in general a twodimensional testcase. We set it up in three dimesnions, because in a common simulation we
will only use three dimensions. So you can check the third dimension, but if your code works proper it should remain zero everywhere */
