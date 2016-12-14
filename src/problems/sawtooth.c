#include "../globals.h"


float Sawtooth_Density(const int ipart)
{
    float x = P[ipart].Pos[0] / Problem.Boxsize[0];
    if (x > 0.5) {
        x -= 0.5;
    }

    double volume = Problem.Boxsize[0]*Problem.Boxsize[1]*Problem.Boxsize[2];
    double mass = Param.Npart * Problem.Mpart;

    double rho = mass/volume;

    return 2.0*rho*x;
}
