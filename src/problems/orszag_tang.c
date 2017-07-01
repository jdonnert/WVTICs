#include "../globals.h"

void setup_Orszag_Tang_Vortex()
{
    Problem.Boxsize[0] = 1.0;
    Problem.Boxsize[1] = 1.0;
    Problem.Boxsize[2] = 0.1;

    sprintf ( Problem.Name, "IC_Orszag_Tang" );

    const double rho = 25 / ( 36 * pi );

    Problem.Rho_Max = rho;
    Problem.Mpart = rho * ( Problem.Boxsize[0] * Problem.Boxsize[1] * Problem.Boxsize[2] ) / Param.Npart;

    Density_Func_Ptr = &Orszag_Tang_Vortex_Density;
    U_Func_Ptr = &Orszag_Tang_Vortex_U;
    Velocity_Func_Ptr = &Orszag_Tang_Vortex_Velocity;
    Magnetic_Field_Func_Ptr = &Orszag_Tang_Vortex_Magnetic_Field;
}

/* At first we set up a constant density in the Box */
float Orszag_Tang_Vortex_Density ( const int ipart )
{
    return 25 / ( 36 * pi );
}

void Orszag_Tang_Vortex_Velocity ( const int ipart, float out[3] )
{

    double const x = Problem.Boxsize[0];
    double const y = Problem.Boxsize[1];
    double const z = Problem.Boxsize[2];

    out[0] = -sin ( 2 * pi * y );
    out[1] = sin ( 2 * pi * x );
    out[2] = 0.0;
}

void Orszag_Tang_Vortex_Magnetic_Field ( const int ipart, float out[3] )
{

    double const x = Problem.Boxsize[0];
    double const y = Problem.Boxsize[1];
    double const z = Problem.Boxsize[2];

    out[0] = -sin ( 2 * pi * y ) / sqrt ( 4 * pi );
    out[1] = sin ( 2 * pi * x ) / sqrt ( 4 * pi );
    out[2] = 0.0;


}

/* We set up the internal energy via the pressure profile and ideal equation of state */

float Orszag_Tang_Vortex_U ( const int ipart )
{
    const double gamma = 5.0 / 3.0;
    const double rho = 25 / ( 36 * pi );
    const double pressure = 5 / 12 * pi;

    return pressure / ( gamma - 1 ) / rho;
}

/* Just a note at the end, the gresho vortex is in general a twodimensional testcase. We set it up in three dimesnions, because in a common simulation we
will only use three dimensions. So you can check the third dimension, but if your code works proper it should remain zero everywhere */

