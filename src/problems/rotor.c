#include "../globals.h"

void setup_Rotor()
{
    Problem.Boxsize[0] = 1.0;
    Problem.Boxsize[1] = 1.0;
    Problem.Boxsize[2] = 0.1;

    sprintf ( Problem.Name, "IC_Rotor" );

    const double rho = 25.0 / ( 36.0 * pi );

    Problem.Rho_Max = rho;
    Problem.Mpart = rho * ( Problem.Boxsize[0] * Problem.Boxsize[1] * Problem.Boxsize[2] ) / Param.Npart;

    Density_Func_Ptr = &Rotor_Density;
    U_Func_Ptr = &Rotor_U;
    Velocity_Func_Ptr = &Rotor_Velocity;
    Magnetic_Field_Func_Ptr = &Rotor_Magnetic_Field;
}

/* At first we set up a constant density in the Box */
float Rotor_Density ( const int ipart )
{
    double Radius = sqrt ( P[ipart].Pos[0] * P[ipart].Pos[0] + P[ipart].Pos[1] * P[ipart].Pos[1] );

    if ( Radius > 0 && Radius <= 0.1 ) {
        return 10.0;
    } else if ( Radius > 0.1 && Radius <= 0.115 ) {
        return 1 + 9 * ( ( 0.115 - Radius ) / ( 0.115 - 0.1 ) );
    } else {
        return 1.0;
    }
}

void Rotor_Velocity ( const int ipart, float out[3] )
{

    double const x = P[ipart].Pos[0] - Problem.Boxsize[0] * 0.5;
    double const y = P[ipart].Pos[1] - Problem.Boxsize[1] * 0.5;

    double Radius = sqrt ( x * x + y * y );

    if ( Radius > 0 && Radius <= 0.1 ) {
        out[0] = -2 * y / 0.1;
        out[1] = 2 * x / 0.1;
        out[2] = 0.0;
    } else if ( Radius > 0.1 && Radius <= 0.115 ) {
        out[0] =  -2 * y / 0.1 * ( ( 0.115 - Radius ) / ( 0.115 - 0.1 ) ) / Radius;
        out[1] =  2 * x / 0.1 * ( ( 0.115 - Radius ) / ( 0.115 - 0.1 ) ) / Radius;
        out[2] = 0.0;
    } else {
        out[0] = 0.0;
        out[1] = 0.0;
        out[2] = 0.0;
    }


}

void Rotor_Magnetic_Field ( const int ipart, float out[3] )
{

    out[0] = 5 / sqrt ( 4 * pi );
    out[1] = 0.0;
    out[2] = 0.0;

}

/* We set up the internal energy via the pressure profile and ideal equation of state */

float Rotor_U ( const int ipart )
{
    const double gamma = 5.0 / 3.0;
    const double rho = 25.0 / ( 36.0 * pi );
    const double pressure = 5.0 / 12.0 * pi;

    return pressure / ( gamma - 1 ) / rho;
}


