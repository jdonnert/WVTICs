#include "../globals.h"

void setup_Sod_Shock()
{
    Problem.Boxsize[0] = 140;
    Problem.Boxsize[1] = 1;
    Problem.Boxsize[2] = 1;

    sprintf ( Problem.Name, "IC_SodShock" );

    Problem.Rho_Max = 1.0;

    Density_Func_Ptr = &Sod_Shock_Density;
    U_Func_Ptr = &Sod_Shock_U;
}

float Sod_Shock_Density ( const int ipart )
{
    const double rhoLeft = 1.0, rhoRight = 0.125;

    if ( P[ipart].Pos[0] <= 0.5 * Problem.Boxsize[0] ) {
        return rhoLeft;
    } else {
        return rhoRight;
    }
}

float Sod_Shock_U ( const int ipart )
{
    const double gamma = 5.0 / 3.0;
    const double rhoLeft = 1.0, rhoRight = 0.125;
    const double pLeft = 1.0, pRight = 0.1;

    if ( P[ipart].Pos[0] <= 0.5 * Problem.Boxsize[0] ) {
        return pLeft / ( gamma - 1.0 ) / rhoLeft;
    } else {
        return pRight / ( gamma - 1.0 ) / rhoRight;
    }
}
