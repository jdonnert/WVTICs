#include "globals.h"

/* Standard Gadget units */
const double ULength = 3.08568025e21; 	// kpc in cgs
const double UMass = 1.989e43;			// 10^10 Msol in cgs
const double UVel = 1e5;					// km/s in cgs

float zero_function(const int ipart);
void zero_function_vec(const int ipart, float out[3]);

void setup_problem(const int Flag, const int Subflag);
void mpart_from_integral();

void Setup()
{
    size_t nBytes = Param.Npart * sizeof(*P); // allocate particles
    P = Malloc(nBytes);
	memset(P, 0, nBytes);

    nBytes = Param.Npart * sizeof(*SphP);
	SphP = Malloc(nBytes);
	memset(SphP, 0, nBytes);

    setup_problem(Param.Problem_Flag, Param.Problem_Subflag);

    mpart_from_integral();

    printf( "Problem %d.%d : %s \n"
                    "   Npart: %d \n"
                    "   Mpart: %g \n"
                    "   Boxsize: %g x %g x %g \n"
                    "   Periodic: %d \n\n",
            Param.Problem_Flag, Param.Problem_Subflag, Problem.Name,
            Param.Npart, Problem.Mpart, Problem.Boxsize[0],
            Problem.Boxsize[1], Problem.Boxsize[2], Problem.Periodic );
}

void setup_problem(const int Flag, const int Subflag)
{
	Density_Func_Ptr = &zero_function;
	U_Func_Ptr = &zero_function;
	Velocity_Func_Ptr = &zero_function_vec;
	Magnetic_Field_Func_Ptr = &zero_function_vec;

	Problem.Mpart = 1; // required to renormalize later

    Problem.Periodic = true; // standard settings
	Problem.Rho_Max = 1.0;
	Problem.Boxsize[0] = Problem.Boxsize[1] = Problem.Boxsize[2] = 1;

	switch (Flag) {

		case 0:

			switch (Subflag) {

				case 0:

                    Problem.Boxsize[0] = 1;
                    Problem.Boxsize[1] = 1;
                    Problem.Boxsize[2] = 0.5;

					sprintf(Problem.Name, "IC_Constant_Density");

					Density_Func_Ptr = &Constant_Density;

					break;

				case 1:
                    
					Problem.Boxsize[0] = 1;
                    Problem.Boxsize[1] = 0.5;
                    Problem.Boxsize[2] = 0.1;

					sprintf(Problem.Name, "IC_TopHat");

					Problem.Rho_Max = 10* Param.Npart * Problem.Mpart 
						/ (Problem.Boxsize[0] * Problem.Boxsize[1] * Problem.Boxsize[2]) ;

					Density_Func_Ptr = &TopHat_Density;

					break;

				case 2:

					sprintf(Problem.Name, "IC_Sawtooth");

					Density_Func_Ptr = &Sawtooth_Density;

					Problem.Rho_Max = 10* Param.Npart * Problem.Mpart 
						/ (Problem.Boxsize[0] * Problem.Boxsize[1] * Problem.Boxsize[2]) ;

					break;

				case 3:

					sprintf(Problem.Name, "IC_SineWave");

					Problem.Boxsize[0] = 1;
                    Problem.Boxsize[1] = 0.75;
                    Problem.Boxsize[2] = 0.75;

					Density_Func_Ptr = &SineWave_Density;

					Problem.Rho_Max = 10* Param.Npart * Problem.Mpart 
						/ (Problem.Boxsize[0] * Problem.Boxsize[1] * Problem.Boxsize[2]) ;

					break;

				default:

					Assert(false,"Effect %d.%d not implemented",Flag, Subflag);

					break;
			}

			break;

        case 1:

            Problem.Periodic = false;

            switch (Subflag) {

                case 0:

                    sprintf(Problem.Name, "IC_GradientDensity");

                    Density_Func_Ptr = &Gradient_Density;

                    break;

                default:

                    Assert(false, "Effect %d.%d not implemented", Flag, Subflag);

                    break;
            }

            break;

		case 2:

			switch (Subflag) {

				case 0:


					sprintf(Problem.Name, "IC_Magneticum");

					Problem.Boxsize[0] = 1;
					Problem.Boxsize[1] = 1;
					Problem.Boxsize[2] = 0.8;

					Density_Func_Ptr = &Magneticum_Density;

					break;


				case 1:
#ifndef EAT_PNG
					printf("Error: must use OPT += -DEAT_PNG in Makefile\n");

					exit(1);
#else

					sprintf(Problem.Name, "IC_PNG");

					Setup_Density_From_Image();

					Density_Func_Ptr = &Png_Density;

					// Problem.Periodic = false;
				        Problem.Boxsize[0] = Image.Xpix;
					Problem.Boxsize[1] = Image.Ypix;
					Problem.Boxsize[2] = Image.Zpix;
#endif // EAT_PNG

					break;

				default:

					Assert(false, "Effect %d.%d not implemented", Flag, Subflag);

					break;
			}

			break;

		case 3:

			sprintf(Problem.Name, "IC_DoubleShock");

			Problem.Boxsize[0] = 2000; // [kpc]
			Problem.Boxsize[1] = 200;
			Problem.Boxsize[2] = 100;

			Problem.Rho_Max = 10* Param.Npart * Problem.Mpart 
					/ (Problem.Boxsize[0] * Problem.Boxsize[1] * Problem.Boxsize[2]) ;

			Density_Func_Ptr = &Double_Shock_Density;
			U_Func_Ptr = &Double_Shock_U;
			Velocity_Func_Ptr = &Double_Shock_Velocity;

			Setup_Double_Shock(Subflag);

			break;

		default:

			Assert(false, "Effect %d.%d not implemented", Flag, Subflag);

			break;
	}

	Assert(Problem.Boxsize[0] >= Problem.Boxsize[1] &&
		   Problem.Boxsize[0] >= Problem.Boxsize[2],
		   "Boxsize[0] has to be largest for ngb finding to work.");

	return ;
}

float zero_function(const int ipart)
{
	return 0;
}

void zero_function_vec(const int ipart, float out[3])
{
	out[0] = out[1] = out[2] = 0;

	return;
}

void mpart_from_integral()
{
    const int N = 1ULL << 9;
    const double dx = Problem.Boxsize[0]/N, 
		  		 dy = Problem.Boxsize[1]/N, 
				 dz = Problem.Boxsize[2]/N;
	
    double tot_mass = 0.0;

    for (int i = 0; i < N; i++)
    {
        P[0].Pos[0] = (i + 0.5) * dx;

        for (int j = 0; j < N; j++)
        {
            P[0].Pos[1] = (j + 0.5) * dy;

            for (int k = 0; k < N; k++)
            {
                P[0].Pos[2] = (k + 0.5) * dz;

                tot_mass += Density_Func_Ptr(0) * dx * dy * dz;
            }
        }
    }

    Problem.Mpart = tot_mass / Param.Npart;

	Assert(Problem.Mpart > 0, 
			"Particle mass has to be finite, have %g", Problem.Mpart);

	return ;
}
