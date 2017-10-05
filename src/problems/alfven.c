#include "../globals.h"

#define ORDER 1
#define RHO_MEAN 1.0

void setup_Linear_Alfven_Wave()
{
    /*
     * Hopkins & Raives 2016 (http://adsabs.harvard.edu/abs/2016MNRAS.455...51H)
     * See Sec. 3.1 (Linear magnetosonic wave), and Fig. 1.
     * Also see Stone+ 2008 (http://adsabs.harvard.edu/abs/2008ApJS..178..137S), Sec. 8.2
     */

    sprintf ( Problem.Name, "IC_LinearAlfvenWave" );

    Problem.Boxsize[0] = 1;
    Problem.Boxsize[1] = 0.1;
    Problem.Boxsize[2] = 0.1;

    Problem.Rho_Max = 1.0 + 1e-6;

    Density_Func_Ptr = &Linear_Alfven_Wave_Density;
}

float Linear_Alfven_Wave_Density ( const int ipart )
{
    float x = P[ipart].Pos[0] / Problem.Boxsize[0];

    return RHO_MEAN * ( 1 + 1e-6 * sin ( 2 * pi * ORDER * x ) );
}


void Linear_Alfven_Wave_Velocity ( const int ipart, float out[3] )
{
    out[0] = out[1] = out[2] = 0.0;

    return;
}

void Linear_Alfven_Wave_Magnetic_Field ( const int ipart, float out[3] )
{
    const double sqrt_fourpi = sqrt ( 4 * pi );

    out[0] = sqrt_fourpi * 1.0;
    out[1] = sqrt_fourpi * sqrt2;
    out[2] = sqrt_fourpi * 0.5;

    return;
}

/* We set up the internal energy via the pressure profile and ideal equation of state */

float Linear_Alfven_Wave_U ( const int ipart )
{
    const double gamma = 5.0 / 3.0;
    double pressure = 1.0 / gamma;

    return pressure / ( gamma - 1 ) / Linear_Alfven_Wave_Density ( ipart );
}
