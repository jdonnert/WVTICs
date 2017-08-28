#include "../globals.h"

void setup_Zeldovich_Pancake()
{
    Problem.Boxsize[0] = 64.0;
    Problem.Boxsize[1] = 64.0;
    Problem.Boxsize[2] = 64.0;

    sprintf ( Problem.Name, "IC_Zeldovich_Pancake" );

    const double Hubble = 67.74; // Planck 2015 value
    const double G =  6.67259e-8;
    const double rho = 3 * Hubble * Hubble / 8 / pi / G; // Critical density of the Universe

    Problem.Rho_Max = rho;

    Problem.Mpart = rho * p3 ( Problem.Boxsize[0] )  / Param.Npart;

    Density_Func_Ptr = &Zeldovich_Pancake_Density;
    U_Func_Ptr = &Zeldovich_Pancake_U;
    Velocity_Func_Ptr = &Zeldovich_Pancake_Velocity;

}


/* At first we set up a constant density in the Box */
float Zeldovich_Pancake_Density ( const int ipart )
{
    const double redshift_start = 100.0;
    const double redshift_crit = 1.0;
    const double Hubble = 67.74;
    const double G = 6.67259e-8;
    const double rho = 3 * Hubble * Hubble / 8 / pi / G;
    const double  k = 2 * pi / Problem.Boxsize[0];

    return rho / ( 1 - ( 1 + redshift_crit ) / ( 1 + redshift_start ) * cos ( k * P[ipart].Pos[0] ) );
}

void Zeldovich_Pancake_Velocity ( const int ipart, float out[3] )
{
    const double redshift_start = 100.0;
    const double redshift_crit = 1.0;
    const double Hubble = 67.74;
    const double k = 2 * pi / Problem.Boxsize[0];

    out[0] = - Hubble * ( 1 + redshift_crit ) / sqrt ( 1 + redshift_start ) * sin ( k * P[ipart].Pos[0] ) / k;
    out[1] = 0.0;
    out[2] = 0.0;
}

/* We set up the internal energy via the pressure profile and ideal equation of state */

float Zeldovich_Pancake_U ( const int ipart )
{
    const double redshift_start = 100.0;
    const double redshift_crit = 1.0;
    const double Hubble = 67.74;
    const double G = 6.67259e-8;
    const double k = 2 * pi / Problem.Boxsize[0];
    const double gamma = 5. / 3.;

    //! @todo set this up correctly
    return 1.0;

}

