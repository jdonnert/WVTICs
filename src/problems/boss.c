#include "../globals.h"

void setup_Boss()
{
    Problem.Boxsize[0] = 5;
    Problem.Boxsize[1] = 5;
    Problem.Boxsize[2] = 5;

    sprintf ( Problem.Name, "IC_Boss" );

    const double rho = 3.82 * 10e-18; // This value is empiric
    const double Radius = 5 * 10e16 ;

    Problem.Rho_Max = rho;

    Problem.Mpart = rho * 4.0 / 3.0 * pi * p3 ( Radius )  / Param.Npart;

    Density_Func_Ptr = &Boss_Density;

}

float Boss_Phi ( const int ipart )
{
    double const x = P[ipart].Pos[0] - Problem.Boxsize[0] * 0.5;
    double const y = P[ipart].Pos[1] - Problem.Boxsize[1] * 0.5;
    double const z = P[ipart].Pos[2] - Problem.Boxsize[2] * 0.5;

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
    const double rho = 3.82 * 10e-18;

    if ( Radius < 1 ) {
        return rho * ( 1 + 0.1 * cos ( 2 * Boss_Phi ( ipart ) ) );
    } else {
        return 0.0;
    }

}

