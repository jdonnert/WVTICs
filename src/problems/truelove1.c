#include "../globals.h"

void setup_Truelove1()
{
    Problem.Periodic[0] = Problem.Periodic[1] = Problem.Periodic[2] = false;

    Problem.Boxsize[0] = 0.032; //5e16 cm in parsec
    Problem.Boxsize[1] = 0.032;
    Problem.Boxsize[2] = 0.032;

    sprintf ( Problem.Name, "IC_Truelove1" );

    const double rho =  234177.1; // 10^-16.8 g/cm^3 in solar masses / pc ^3 

    Problem.Rho_Max = rho * 1.1;

    Density_Func_Ptr = &Truelove1_Density;

}

float Truelove1_Phi ( double const x, double const y )
{

   return atan2( y , x );

}

/* At first we set up a constant density in the Box */
float Truelove1_Density ( const int ipart , const double bias )
{
    double const x = P[ipart].Pos[0] - Problem.Boxsize[0] * 0.5;
    double const y = P[ipart].Pos[1] - Problem.Boxsize[1] * 0.5;
    double const z = P[ipart].Pos[2] - Problem.Boxsize[2] * 0.5;

    double Radius = sqrt ( x*x + y*y + z*z );
    double R1 = 9.398e-3;
    const double rho = 234177.1*exp(-pow(Radius/R1, 2));

    return rho * ( 1 + 0.1 * cos ( 2 * Truelove1_Phi ( x, y ) ) );

}
