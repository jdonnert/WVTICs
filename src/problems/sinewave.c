#include "../globals.h"

#define ORDER 4

float SineWave_Density(const int ipart)
{
    float x = P[ipart].Pos[0] / Problem.Boxsize[0];

    double volume = Problem.Boxsize[0]*Problem.Boxsize[1]*Problem.Boxsize[2];
    double mass = Param.Npart * Problem.Mpart;

    const double rho0 = mass/volume;

    return rho0 * (1 + 0.5 * sin(2 * pi * ORDER * x));
}
