#include "../globals.h"

void setup_Strong_Blast()
{
    /*
     * Hopkins & Raives 2016 (http://adsabs.harvard.edu/abs/2016MNRAS.455...51H)
     * See Sec. 3.8, and Fig. 17 & 18.
     */

    Problem.Boxsize[0] = 1.0;
    Problem.Boxsize[1] = 1.0;
    Problem.Boxsize[2] = 0.1;

    sprintf ( Problem.Name, "IC_StrongBlast" );

    const double rho = 1.0;

    Problem.Rho_Max = rho;

    Density_Func_Ptr = &Strong_Blast_Density;
    U_Func_Ptr = &Strong_Blast_U;
    Velocity_Func_Ptr = &Strong_Blast_Velocity;
    Magnetic_Field_Func_Ptr = &Strong_Blast_Magnetic_Field;
}

/* At first we set up a constant density in the Box */
float Strong_Blast_Density ( const int ipart )
{
    return 1.0;
}

void Strong_Blast_Velocity ( const int ipart, float out[3] )
{
    out[0] = out[1] = out[2] = 0.0;

    return;
}

void Strong_Blast_Magnetic_Field ( const int ipart, float out[3] )
{

    out[0] = 1 / sqrt ( 2 );
    out[1] = 1 / sqrt ( 2 );
    out[2] = 0.0;

    return;
}

/* We set up the internal energy via the pressure profile and ideal equation of state */

float Strong_Blast_U ( const int ipart )
{

    double const x = P[ipart].Pos[0] - Problem.Boxsize[0] * 0.5;
    double const y = P[ipart].Pos[1] - Problem.Boxsize[1] * 0.5;

    double radius = sqrt ( p2 ( x ) + p2 ( y ) );

    const double gamma = 5.0 / 3.0;
    double pressure = 1.0;


    if ( radius > 0 && radius <= 0.1 ) {
        pressure = 10.0;
    } else if ( radius > 0.1 ) {
        pressure = 0.1;
    }


    return pressure / ( gamma - 1 ) / Strong_Blast_Density ( ipart );
}
