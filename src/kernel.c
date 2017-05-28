#include "kernel.h"
#include "globals.h"

float sph_kernel ( const float r, const float h )
{
#ifdef SPH_CUBIC_SPLINE
    return sph_kernel_M4 ( r, h );
#else
#ifdef SPH_WC2
    return sph_kernel_WC2 ( r, h );
#else
    return sph_kernel_WC6 ( r, h );
#endif // SPH_WC2
#endif // SPH_CUBIC_SPLINE
}

float sph_kernel_derivative ( const float r, const float h )
{
#ifdef SPH_CUBIC_SPLINE
    return sph_kernel_derivative_M4 ( r, h );
#else
#ifdef SPH_WC2
    return sph_kernel_derivative_WC2 ( r, h );
#else
    return sph_kernel_derivative_WC6 ( r, h );
#endif // SPH_WC2
#endif // SPH_CUBIC_SPLINE
}

double bias_correction ( const float h )
{
#ifndef SPH_CUBIC_SPLINE
#ifdef SPH_WC2
    return bias_correction_WC2 ( h );
#else
    return bias_correction_WC6 ( h );
#endif // SPH_WC2
#else
    return 0.0;
#endif  // SPH_CUBIC_SPLINE
}

float sph_kernel_WC6 ( const float r, const float h )
{
    const double u = r / h;
    const double t = 1 - u;

    return 1365.0 / ( 64 * pi ) / p3 ( h ) * t * t * t * t * t * t * t * t * ( 1 + 8 * u + 25 * u * u + 32 * u * u * u );
}

float sph_kernel_derivative_WC6 ( const float r, const float h )
{
    const float u = r / h;
    const double t = 1 - u;

    return 1365.0 / ( 64 * pi ) / ( h * h * h * h ) * -22.0 * t * t * t * t * t * t * t * u * ( 16 * u * u + 7 * u + 1 );
}

double bias_correction_WC6 ( const float h )
{
    return -0.0116 * pow ( DESNNGB * 0.01, -2.236 ) * Problem.Mpart * sph_kernel_WC6 ( 0, h );
}

float sph_kernel_WC2 ( const float r, const float h )
{
    const double u = r / h;
    const double t = 1 - u;

    return 21.0 / ( 2 * pi ) / p3 ( h ) * t * t * t * t * ( 1 + 4 * u );
}

float sph_kernel_derivative_WC2 ( const float r, const float h )
{
    const float u = r / h;
    const double t = 1 - u;

    return 21.0 / ( 2 * pi ) / ( h * h * h * h ) * -20.0 * t * t * t * u;
}

double bias_correction_WC2 ( const float h )
{
    return -0.0294 * pow ( DESNNGB * 0.01, -0.977 ) * Problem.Mpart * sph_kernel_WC2 ( 0, h ); // WC2 (Dehnen+ 12)
}

float sph_kernel_M4 ( const float r, const float h ) // cubic spline
{
    double wk = 0;
    double u = r / h;

    if ( u < 0.5 ) {
        wk = ( 2.546479089470 + 15.278874536822 * ( u - 1 ) * u * u );
    } else {
        wk = 5.092958178941 * ( 1.0 - u ) * ( 1.0 - u ) * ( 1.0 - u );
    }

    return wk / p3 ( h );
}

float sph_kernel_derivative_M4 ( const float r, const float h )
{
    double dwk = 0;
    double u = r / h;

    if ( u < 0.5 ) {
        dwk = u * ( 45.836623610466 * u - 30.557749073644 );
    } else {
        dwk = ( -15.278874536822 ) * ( 1.0 - u ) * ( 1.0 - u );
    }

    return dwk / ( h * h * h * h );
}
