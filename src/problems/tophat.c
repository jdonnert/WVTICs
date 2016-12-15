#include "../globals.h"

#define DENSITY_STEP 0.5

float TopHat_Density(const int ipart)
{
    float x = P[ipart].Pos[0] / Problem.Boxsize[0];

    double volume = Problem.Boxsize[0]*Problem.Boxsize[1]*Problem.Boxsize[2];
    double mass = Param.Npart * Problem.Mpart;

    const double rho0 = mass/volume;
    const double halfstep = DENSITY_STEP * rho0;

    double rho_max = rho0 + halfstep;
    double rho_min = rho0 - halfstep;

    if (x < 0.1 || x > 0.9) {
        return rho_min;
    }
    if (x > 0.4 && x < 0.6) {
        return rho_max;
    }
    if (x > 0.6) {
        return rho_max - (rho_max - rho_min) * (x - 0.6) / 0.3;
    }
    return rho_min + (rho_max - rho_min) * (x - 0.4) / 0.3;
}
