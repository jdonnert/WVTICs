#include "globals.h"

void Make_Positions()
{
    printf("Sampling positions ..."); fflush(stdout);

	#pragma omp parallel for
    for (int ipart = 0; ipart < Param.Npart; ipart++) {

#ifdef REJECTION_SAMPLING
        double rho = 0.0, rho_r = 0.0;

        while (rho >= rho_r) {

            P[ipart].Pos[0] = erand48(Omp.Seed) * Problem.Boxsize[0];
            P[ipart].Pos[1] = erand48(Omp.Seed) * Problem.Boxsize[1];
            P[ipart].Pos[2] = erand48(Omp.Seed) * Problem.Boxsize[2];

			rho = Problem.Rho_Max * erand48(Omp.Seed);
			rho_r = Density_Func_Ptr(ipart);
        }
#else
        P[ipart].Pos[0] = erand48(Omp.Seed) * Problem.Boxsize[0];
        P[ipart].Pos[1] = erand48(Omp.Seed) * Problem.Boxsize[1];
        P[ipart].Pos[2] = erand48(Omp.Seed) * Problem.Boxsize[2];
#endif //REJECTION_SAMPLING

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
		SphP[ipart].Bfld[0] = 0.0;
		SphP[ipart].Bfld[1] = 0.0;
		SphP[ipart].Bfld[2] = 0.0;
    }

    printf(" done\n");

	return ;
}

