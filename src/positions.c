#include "globals.h"

void Make_Positions()
{
    printf("Sampling positions ..."); fflush(stdout);
	
	#pragma omp parallel for
    for (int ipart = 0; ipart < Param.Npart; ipart++) {

        for (;;) {

            P[ipart].Pos[0] = erand48(Omp.Seed);
            P[ipart].Pos[1] = erand48(Omp.Seed);
            P[ipart].Pos[2] = erand48(Omp.Seed);

			double rho = Problem.Rho_Max * erand48(Omp.Seed);

			double rho_r = Density_Func_Ptr(ipart);

			if (rho < rho_r) // rejection sampling
				break;
        }
		
		P[ipart].Type = 0;
    }

    printf(" done\n");

    return;
}

void Make_Velocities()
{

    printf("Velocities ..."); fflush(stdout);

	#pragma omp parallel for
    for (int ipart = 0; ipart < Param.Npart; ipart++) {
	
		(*Velocity_Func_Ptr) (ipart, P[ipart].Vel);
	}

    printf(" done\n");

	return ;
}

void Make_Temperatures()
{

    printf("Internal Energy ..."); fflush(stdout);

	#pragma omp parallel for
    for (int ipart = 0; ipart < Param.Npart; ipart++) {

		SphP[ipart].U = (*U_Func_Ptr) (ipart);
	}

    printf(" done\n");

	return ;
}

void Make_Magnetic_Fields()
{

    printf("Magnetic Field ..."); fflush(stdout);

	#pragma omp parallel for
    for (int ipart = 0; ipart < Param.Npart; ipart++) {
	
		(*Magnetic_Field_Func_Ptr) (ipart, SphP[ipart].Bfld);
	}

    printf(" done\n");

	return ;
}

