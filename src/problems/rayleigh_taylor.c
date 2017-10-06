#include "../globals.h"

void setup_Rayleigh_Taylor_Instability()
{
    /*
     * Hopkins & Raives 2016 (http://adsabs.harvard.edu/abs/2016MNRAS.455...51H)
     * See Sec. 3.9, and Fig. 19 & 20
     *
     * CAUTION, we changed x and y with respect to the paper because Boxsize[0]
     * has to be largest for ngb finding to work.
     * CAUTION2, the paper specifies x has periodic boundaries; y reflecting.
     * Here we have periodic boundaries for both x and y (but perhaps we require
     * reflective boundaries in the x-direction (because x/y-swap)
     *
     * http://www.tapir.caltech.edu/~phopkins/sims/rt.params
     *
     * Also see Abel+ 2011 (http://adsabs.harvard.edu/abs/2011MNRAS.413..271A),
     */

    sprintf ( Problem.Name, "IC_RayleighTaylorInstability" );

    Problem.Boxsize[0] = 1.0;
    Problem.Boxsize[1] = 0.5;
    Problem.Boxsize[2] = 0.1;

    Problem.Rho_Max = 2.0;

    Density_Func_Ptr = &Rayleigh_Taylor_Instability_Density;
    U_Func_Ptr = &Rayleigh_Taylor_Instability_U;
    Velocity_Func_Ptr = &Rayleigh_Taylor_Instability_Velocity;
    Magnetic_Field_Func_Ptr = &Rayleigh_Taylor_Instability_Magnetic_Field;
}

float Rayleigh_Taylor_Instability_Density ( const int ipart )
{
    float x = P[ipart].Pos[0] / Problem.Boxsize[0];
    const float rho1 = 1.0f;
    const float rho2 = 2.0f;
    const float delta = 0.025f;

    return rho1 + ( rho2 - rho1 ) / ( 1 + exp ( - ( x - 0.5 ) / delta ) );
}


void Rayleigh_Taylor_Instability_Velocity ( const int ipart, float out[3] )
{
    float x = P[ipart].Pos[0] / Problem.Boxsize[0];
    float y = P[ipart].Pos[1] / Problem.Boxsize[1];

    if ( x > 0.3 && x < 0.7 ) {
        // Density perturbation in the x-direction
        out[0] = 0.025 * ( 1 + cos ( 8 * pi * ( y + 0.25 ) ) ) * ( 1 + sin ( 5 * pi * ( x - 0.5 ) ) );
    } else {
        out[0] = 0.0;
    }

    out[1] = out[2] = 0.0;



    return;
}

void Rayleigh_Taylor_Instability_Magnetic_Field ( const int ipart, float out[3] )
{
    const float sqrt_fourpi = sqrt ( 4 * pi );

    out[0] = 0.0;
    out[1] = sqrt_fourpi * 0.07;
    out[2] = 0.0;

    return;
}

/* We set up the internal energy via the pressure profile and ideal equation of state */

float Rayleigh_Taylor_Instability_U ( const int ipart )
{
    const float gamma = 1.4;
    const float rho1 = 1.0;
    const float rho2 = 2.0;
    const float delta = 0.025;
    const float grav_acc = -0.5;

    /* "Initial pressure is assigned to produce Hydrostatic Equilibrium with a
     * uniform gravitational acceleration g = -1/2 in the y-direction (at the
     * interface, P = rho2/gamma = 10/7 so the soundspeed c_s = 1" - taken from
     * Hopkins & Raives 2016 Sec. 3.9, second paragraph
     *
     * Hydrostatic Equation: rho F = nabla P --> dP/dx = rho g --> P = P_0 + rho g x
     * (x because we swapped y and x here)
     */

    float x = P[ipart].Pos[0] / Problem.Boxsize[0];
    float pressure = rho2 / gamma + grav_acc * Rayleigh_Taylor_Instability_Density ( ipart ) * ( x - 0.5 );

    return pressure / ( gamma - 1 ) / Rayleigh_Taylor_Instability_Density ( ipart );
}
