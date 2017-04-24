#include "../globals.h"

#define ORDER 4

float SineWave_Density(const int ipart)
{
    float x = P[ipart].Pos[0] / Problem.Boxsize[0];

    const double rho0 = 1.0;

    return rho0 * (1 + 0.5 * sin(2 * pi * ORDER * x));
}
