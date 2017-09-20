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

    Density_Func_Ptr = &Zeldovich_Pancake_Density;
    U_Func_Ptr = &Zeldovich_Pancake_U;
    Velocity_Func_Ptr = &Zeldovich_Pancake_Velocity;

}

float q_of_x ( const int ipart )
{

    return 0.0120544 + 0.999977 * P[ipart].Pos[0];
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

    return rho / ( 1 - ( 1 + redshift_crit ) / ( 1 + redshift_start ) * cos ( k * q_of_x ( ipart ) ) );
}

void Zeldovich_Pancake_Velocity ( const int ipart, float out[3] )
{
    const double redshift_start = 100.0;
    const double redshift_crit = 1.0;
    const double Hubble = 67.74;
    const double k = 2 * pi / Problem.Boxsize[0];

    out[0] = - Hubble * ( 1 + redshift_crit ) / sqrt ( 1 + redshift_start ) * sin ( k * q_of_x ( ipart ) ) / k;
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
    const double gamma = 5. / 3.;
    const double rho = 3 * Hubble * Hubble / 8 / pi / G;
    const double temp_zero = 1.0;
    const double kboltzmann = 1.380658e-16;
    const double yhelium = ( 1 - 0.76 ) / ( 4 * 0.76 );
    const double mean_mol_weight = ( 1 + 4 * yhelium ) / ( 1 + 3 * yhelium + 1 );
    const double v_unit = 1e5;
    const double prtn = 1.672623e-24;
    const double u_fac = kboltzmann  / ( ( gamma - 1 ) * v_unit * v_unit * prtn * mean_mol_weight );

    return u_fac * temp_zero * pow ( redshift_start / redshift_crit, 2 ) * pow ( SphP[ipart].Rho / rho, 2. / 3. );

}

