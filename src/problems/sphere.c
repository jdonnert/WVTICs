#include "../globals.h"

void setup_Sphere()
{
    Problem.Boxsize[0] = 1;
    Problem.Boxsize[1] = 1;
    Problem.Boxsize[2] = 1;

    sprintf ( Problem.Name, "IC_Sphere" );

    Problem.Rho_Max = 1.0;

    Density_Func_Ptr = &Sphere_Density;
}

float Sphere_Density ( const int ipart , const double bias )
{
    const double rhoIn = 1.0, rhoOut = 0.1;

    const double radius = 0.5;

    double rx,ry,rz;

    rx = P[ipart].Pos[0] - 0.5 * Problem.Boxsize[0];
    ry = P[ipart].Pos[1] - 0.5 * Problem.Boxsize[1];
    rz = P[ipart].Pos[2] - 0.5 * Problem.Boxsize[2];

    double r = sqrt(rx*rx + ry*ry + rz*rz);

    if ( r <= radius) {
        return rhoIn;
    }
    else {
        return rhoOut;
    }
}
