#include "../globals.h"

void setup_Boss()
{
    Problem.Periodic = 0;

    Problem.Boxsize[0] = 0.016; //5e16 cm in parsec
    Problem.Boxsize[1] = 0.016;
    Problem.Boxsize[2] = 0.016;

    sprintf ( Problem.Name, "IC_Boss" );

    const double rho = 5624.3; // 3.82 *10e-18 in g/cm^3 in solar masses per cubic parsec

    Problem.Rho_Max = rho * 1.1;

    Density_Func_Ptr = &Boss_Density;

}

float Boss_Phi ( double const x, double const y )
{
    if ( x > 0 ) {
        return atan ( y / x );
    } else if ( x == 0 ) {
        return y / abs ( y ) * 0.5 * pi;
    } else if ( x < 0 && y >= 0 ) {
        return atan ( y / x ) + pi;
    } else {
        return atan ( y / x ) - pi;
    }
}

/* At first we set up a constant density in the Box */
float Boss_Density ( const int ipart )
{
    double const x = P[ipart].Pos[0] - Problem.Boxsize[0] * 0.5;
    double const y = P[ipart].Pos[1] - Problem.Boxsize[1] * 0.5;
    double const z = P[ipart].Pos[2] - Problem.Boxsize[2] * 0.5;

    double Radius = sqrt ( x * x + y * y + z * z );
    const double rho = 5624.3;

    if ( Radius < 0.016 * 0.5 ) {
        return rho * ( 1 + 0.1 * cos ( 2 * Boss_Phi ( x, y ) ) );
    } else {
        return rho * 0.01;
    }

}

