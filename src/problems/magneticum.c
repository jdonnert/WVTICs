#include "../globals.h"

float Magneticum_Density(const int ipart)
{
	const float x = P[ipart].Pos[0] / Problem.Boxsize[0];
	const float y = P[ipart].Pos[1] / Problem.Boxsize[1];
	const float z = P[ipart].Pos[2] / Problem.Boxsize[2];

	double rho = 1;

    // Square box in the xy plane
    // if ( (4 < x && x < 6) &&
    //      (1 < y && y < 9) &&
    //      (2 < z && z < 8) )
    //     return rho;
    // else
    //     return rho_bkg;
    // return rho_bkg;

    // Magneticum :-) ...
    // M
    if ((0.00*boxsize < x && x < 0.02*boxsize) &&
        (0.45*boxsize < y && y < 0.56*boxsize) &&
        (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    else if ((0.04*boxsize < x && x < 0.06*boxsize) &&
          (0.45*boxsize < y && y < 0.55*boxsize) &&
          (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    else if ((0.08*boxsize < x && x < 0.10*boxsize) &&
          (0.45*boxsize < y && y < 0.55*boxsize) &&
          (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    else if ((0.00*boxsize < x && x < 0.10*boxsize) &&
          (0.53*boxsize < y && y < 0.55*boxsize) &&
          (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    // A
    if ((0.12*boxsize < x && x < 0.14*boxsize) &&
        (0.45*boxsize < y && y < 0.56*boxsize) &&
        (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    else if ((0.17*boxsize < x && x < 0.19*boxsize) &&
          (0.45*boxsize < y && y < 0.55*boxsize) &&
          (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    else if ((0.12*boxsize < x && x < 0.19*boxsize) &&
          (0.53*boxsize < y && y < 0.55*boxsize) &&
          (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    else if ((0.12*boxsize < x && x < 0.19*boxsize) &&
          (0.49*boxsize < y && y < 0.51*boxsize) &&
          (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    // G
    if ((0.21*boxsize < x && x < 0.23*boxsize) &&
        (0.45*boxsize < y && y < 0.55*boxsize) &&
        (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    else if ((0.21*boxsize < x && x < 0.28*boxsize) &&
          (0.53*boxsize < y && y < 0.55*boxsize) &&
          (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    else if ((0.21*boxsize < x && x < 0.28*boxsize) &&
          (0.45*boxsize < y && y < 0.47*boxsize) &&
          (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    else if ((0.26*boxsize < x && x < 0.28*boxsize) &&
          (0.45*boxsize < y && y < 0.49*boxsize) &&
          (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    // N
    if ((0.30*boxsize < x && x < 0.32*boxsize) &&
        (0.45*boxsize < y && y < 0.56*boxsize) &&
        (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    else if ((0.36*boxsize < x && x < 0.38*boxsize) &&
          (0.45*boxsize < y && y < 0.55*boxsize) &&
          (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    else if ((0.30*boxsize < x && x < 0.38*boxsize) &&
         (0.53*boxsize < y && y < 0.55*boxsize) &&
         (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    // E
    if ((0.40*boxsize < x && x < 0.42*boxsize) &&
        (0.45*boxsize < y && y < 0.55*boxsize) &&
        (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    else if ((0.40*boxsize < x && x < 0.48*boxsize) &&
          (0.53*boxsize < y && y < 0.55*boxsize) &&
          (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    else if ((0.40*boxsize < x && x < 0.48*boxsize) &&
          (0.49*boxsize < y && y < 0.51*boxsize) &&
          (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    else if ((0.40*boxsize < x && x < 0.48*boxsize) &&
          (0.45*boxsize < y && y < 0.47*boxsize) &&
          (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    // T
    if ((0.53*boxsize < x && x < 0.55*boxsize) &&
        (0.45*boxsize < y && y < 0.55*boxsize) &&
        (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    else if ((0.50*boxsize < x && x < 0.58*boxsize) &&
          (0.53*boxsize < y && y < 0.55*boxsize) &&
          (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    // I
    if ((0.63*boxsize < x && x < 0.65*boxsize) &&
        (0.45*boxsize < y && y < 0.55*boxsize) &&
        (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    else if ((0.60*boxsize < x && x < 0.68*boxsize) &&
          (0.53*boxsize < y && y < 0.55*boxsize) &&
          (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    else if ((0.60*boxsize < x && x < 0.68*boxsize) &&
          (0.45*boxsize < y && y < 0.47*boxsize) &&
          (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    // C
    if ((0.70*boxsize < x && x < 0.72*boxsize) &&
        (0.45*boxsize < y && y < 0.55*boxsize) &&
        (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    else if ((0.70*boxsize < x && x < 0.78*boxsize) &&
          (0.53*boxsize < y && y < 0.55*boxsize) &&
          (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    else if ((0.70*boxsize < x && x < 0.78*boxsize) &&
          (0.45*boxsize < y && y < 0.47*boxsize) &&
          (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    // U
    if ((0.80*boxsize < x && x < 0.82*boxsize) &&
        (0.45*boxsize < y && y < 0.55*boxsize) &&
        (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    if ((0.86*boxsize < x && x < 0.88*boxsize) &&
        (0.45*boxsize < y && y < 0.55*boxsize) &&
        (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    if ((0.80*boxsize < x && x < 0.88*boxsize) &&
        (0.45*boxsize < y && y < 0.47*boxsize) &&
        (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    // M
    if ((0.90*boxsize < x && x < 0.92*boxsize) &&
        (0.45*boxsize < y && y < 0.56*boxsize) &&
        (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    else if ((0.94*boxsize < x && x < 0.96*boxsize) &&
          (0.45*boxsize < y && y < 0.55*boxsize) &&
          (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    else if ((0.98*boxsize < x && x < 1.00*boxsize) &&
          (0.45*boxsize < y && y < 0.55*boxsize) &&
          (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    else if ((0.90*boxsize < x && x < 1.00*boxsize) &&
          (0.53*boxsize < y && y < 0.55*boxsize) &&
          (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    // Underline
    if ((0.00*boxsize < x && x < 1.00*boxsize) &&
        (0.41*boxsize < y && y < 0.43*boxsize) &&
        (0.0*boxsize < z && z < 1.0*boxsize))
        return rho;
    
	
    return 0.1; // Default
}
