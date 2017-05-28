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

// See [Dehnen & Ali 12]

float sph_kernel_WC6 ( const float r, const float h )
{
    const double u = r / h;
    const double t = 1 - u;

#ifdef TWO_DIM
    double norm = 78.0 / ( 7 * pi ) / p2 ( h );
#else
    double norm = 1365.0 / ( 64 * pi ) / p3 ( h );
#endif //TWO_DIM

    return norm * t * t * t * t * t * t * t * t * ( 1 + 8 * u + 25 * u * u + 32 * u * u * u );
}

float sph_kernel_derivative_WC6 ( const float r, const float h )
{
    const float u = r / h;
    const double t = 1 - u;

#ifdef TWO_DIM
    double norm = 78.0 / ( 7 * pi ) / p2 ( h );
#else
    double norm = 1365.0 / ( 64 * pi ) / p3 ( h );
#endif //TWO_DIM

    return norm / h * -22.0 * t * t * t * t * t * t * t * u * ( 16 * u * u + 7 * u + 1 );
}

double bias_correction_WC6 ( const float h )
{
    return -0.0116 * pow ( DESNNGB * 0.01, -2.236 ) * Problem.Mpart * sph_kernel_WC6 ( 0, h );
}

float sph_kernel_WC2 ( const float r, const float h )
{
    const double u = r / h;
    const double t = 1 - u;

#ifdef TWO_DIM
    double norm = 7.0 / pi / p2 ( h );
#else
    double norm = 21.0 / ( 2 * pi ) / p3 ( h );
#endif //TWO_DIM

    return norm * t * t * t * t * ( 1 + 4 * u );
}

float sph_kernel_derivative_WC2 ( const float r, const float h )
{
    const float u = r / h;
    const double t = 1 - u;

#ifdef TWO_DIM
    double norm = 7.0 / pi / p2 ( h );
#else
    double norm = 21.0 / ( 2 * pi ) / p3 ( h );
#endif //TWO_DIM

    return norm / h * -20.0 * t * t * t * u;
}

double bias_correction_WC2 ( const float h )
{
    return -0.0294 * pow ( DESNNGB * 0.01, -0.977 ) * Problem.Mpart * sph_kernel_WC2 ( 0, h );
}

float sph_kernel_M4 ( const float r, const float h ) // cubic spline
{
    const double u = r / h;
    const double t = 1 - u;
    const double v = ( u > 0.5 ? 0.0 : 0.5 - u );

#ifdef TWO_DIM
    double norm = 80.0 / ( 7.0 * pi ) / p2 ( h );
#else
    double norm = 16.0 / pi / p3 ( h );
#endif //TWO_DIM

    return norm * ( p3 ( t ) - 4.0 * p3 ( v ) );
}

float sph_kernel_derivative_M4 ( const float r, const float h )
{
    const double u = r / h;
    const double t = 1 - u;
    const double v = ( u > 0.5 ? 0.0 : 0.5 - u );

#ifdef TWO_DIM
    double norm = 80.0 / ( 7.0 * pi ) / p2 ( h );
#else
    double norm = 16.0 / pi / p3 ( h );
#endif //TWO_DIM

    return norm / h * ( -3.0 * p2 ( t ) + 12.0 * p2 ( v ) );
}
