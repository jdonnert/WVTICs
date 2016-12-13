#include "../globals.h"

#define DENSITY_CONTRAST 10

float Magneticum_Density(const int ipart)
{
	const float x = P[ipart].Pos[0] / Problem.Boxsize[0];
	const float y = P[ipart].Pos[1] / Problem.Boxsize[1];
	const float z = P[ipart].Pos[2] / Problem.Boxsize[2];
	
	double volume = Problem.Boxsize[0]*Problem.Boxsize[1]*Problem.Boxsize[2];
	double mass = Param.Npart * Problem.Mpart;

	double rho = mass/volume;

    if ((0.00 < x && x < 0.02) &&
        (0.45 < y && y < 0.56) &&
        (0.0 < z && z < 1.0))
        return rho;
    else if ((0.04 < x && x < 0.06) &&
          (0.45 < y && y < 0.55) &&
          (0.0 < z && z < 1.0))
        return rho;
    else if ((0.08 < x && x < 0.10) &&
          (0.45 < y && y < 0.55) &&
          (0.0 < z && z < 1.0))
        return rho;
    else if ((0.00 < x && x < 0.10) &&
          (0.53 < y && y < 0.55) &&
          (0.0 < z && z < 1.0))
        return rho;
    // A
    if ((0.12 < x && x < 0.14) &&
        (0.45 < y && y < 0.56) &&
        (0.0 < z && z < 1.0))
        return rho;
    else if ((0.17 < x && x < 0.19) &&
          (0.45 < y && y < 0.55) &&
          (0.0 < z && z < 1.0))
        return rho;
    else if ((0.12 < x && x < 0.19) &&
          (0.53 < y && y < 0.55) &&
          (0.0 < z && z < 1.0))
        return rho;
    else if ((0.12 < x && x < 0.19) &&
          (0.49 < y && y < 0.51) &&
          (0.0 < z && z < 1.0))
        return rho;
    // G
    if ((0.21 < x && x < 0.23) &&
        (0.45 < y && y < 0.55) &&
        (0.0 < z && z < 1.0))
        return rho;
    else if ((0.21 < x && x < 0.28) &&
          (0.53 < y && y < 0.55) &&
          (0.0 < z && z < 1.0))
        return rho;
    else if ((0.21 < x && x < 0.28) &&
          (0.45 < y && y < 0.47) &&
          (0.0 < z && z < 1.0))
        return rho;
    else if ((0.26 < x && x < 0.28) &&
          (0.45 < y && y < 0.49) &&
          (0.0 < z && z < 1.0))
        return rho;
    // N
    if ((0.30 < x && x < 0.32) &&
        (0.45 < y && y < 0.56) &&
        (0.0 < z && z < 1.0))
        return rho;
    else if ((0.36 < x && x < 0.38) &&
          (0.45 < y && y < 0.55) &&
          (0.0 < z && z < 1.0))
        return rho;
    else if ((0.30 < x && x < 0.38) &&
         (0.53 < y && y < 0.55) &&
         (0.0 < z && z < 1.0))
        return rho;
    // E
    if ((0.40 < x && x < 0.42) &&
        (0.45 < y && y < 0.55) &&
        (0.0 < z && z < 1.0))
        return rho;
    else if ((0.40 < x && x < 0.48) &&
          (0.53 < y && y < 0.55) &&
          (0.0 < z && z < 1.0))
        return rho;
    else if ((0.40 < x && x < 0.48) &&
          (0.49 < y && y < 0.51) &&
          (0.0 < z && z < 1.0))
        return rho;
    else if ((0.40 < x && x < 0.48) &&
          (0.45 < y && y < 0.47) &&
          (0.0 < z && z < 1.0))
        return rho;
    // T
    if ((0.53 < x && x < 0.55) &&
        (0.45 < y && y < 0.55) &&
        (0.0 < z && z < 1.0))
        return rho;
    else if ((0.50 < x && x < 0.58) &&
          (0.53 < y && y < 0.55) &&
          (0.0 < z && z < 1.0))
        return rho;
    // I
    if ((0.63 < x && x < 0.65) &&
        (0.45 < y && y < 0.55) &&
        (0.0 < z && z < 1.0))
        return rho;
    else if ((0.60 < x && x < 0.68) &&
          (0.53 < y && y < 0.55) &&
          (0.0 < z && z < 1.0))
        return rho;
    else if ((0.60 < x && x < 0.68) &&
          (0.45 < y && y < 0.47) &&
          (0.0 < z && z < 1.0))
        return rho;
    // C
    if ((0.70 < x && x < 0.72) &&
        (0.45 < y && y < 0.55) &&
        (0.0 < z && z < 1.0))
        return rho;
    else if ((0.70 < x && x < 0.78) &&
          (0.53 < y && y < 0.55) &&
          (0.0 < z && z < 1.0))
        return rho;
    else if ((0.70 < x && x < 0.78) &&
          (0.45 < y && y < 0.47) &&
          (0.0 < z && z < 1.0))
        return rho;
    // U
    if ((0.80 < x && x < 0.82) &&
        (0.45 < y && y < 0.55) &&
        (0.0 < z && z < 1.0))
        return rho;
    if ((0.86 < x && x < 0.88) &&
        (0.45 < y && y < 0.55) &&
        (0.0 < z && z < 1.0))
        return rho;
    if ((0.80 < x && x < 0.88) &&
        (0.45 < y && y < 0.47) &&
        (0.0 < z && z < 1.0))
        return rho;
    // M
    if ((0.90 < x && x < 0.92) &&
        (0.45 < y && y < 0.56) &&
        (0.0 < z && z < 1.0))
        return rho;
    else if ((0.94 < x && x < 0.96) &&
          (0.45 < y && y < 0.55) &&
          (0.0 < z && z < 1.0))
        return rho;
    else if ((0.98 < x && x < 1.00) &&
          (0.45 < y && y < 0.55) &&
          (0.0 < z && z < 1.0))
        return rho;
    else if ((0.90 < x && x < 1.00) &&
          (0.53 < y && y < 0.55) &&
          (0.0 < z && z < 1.0))
        return rho;
    // Underline
    if ((0.00 < x && x < 1.00) &&
        (0.41 < y && y < 0.43) &&
        (0.0 < z && z < 1.0))
        return rho;
    
	
    return rho / DENSITY_CONTRAST; // Default
}
