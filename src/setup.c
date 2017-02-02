#include "globals.h"

float zero_function(const int ipart);
void zero_function_vec(const int ipart, float out[3]);

void setup_problem(const int Flag, const int Subflag);
void calculateParticleMassByIntegration();

void Setup()
{
    size_t nBytes = Param.Npart * sizeof(*P); // allocate particles
    P = Malloc(nBytes);
	memset(P, 0, nBytes);

    nBytes = Param.Npart * sizeof(*SphP);
	SphP = Malloc(nBytes);
	memset(SphP, 0, nBytes);

    setup_problem(Param.Problem_Flag, Param.Problem_Subflag);
    calculateParticleMassByIntegration();

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

    Problem.Periodic = true; // standard settings
	Problem.Rho_Max = 1.0;
	Problem.Boxsize[0] = Problem.Boxsize[1] = Problem.Boxsize[2] = 1;
	Problem.Mpart = 1.0 / Param.Npart;

	switch (Flag) {

		case 0:

			switch (Subflag) {

				case 0:

                    Problem.Boxsize[0] = 1;
                    Problem.Boxsize[1] = 1;
                    Problem.Boxsize[2] = 0.5;

                    Problem.Mpart = 1.0 / Param.Npart * 
						(Problem.Boxsize[0] * Problem.Boxsize[1] * Problem.Boxsize[2]);

					sprintf(Problem.Name, "IC_Constant_Density");

					Density_Func_Ptr = &Constant_Density;

					break;

				case 1:

					sprintf(Problem.Name, "IC_TopHat");

					Density_Func_Ptr = &TopHat_Density;

					break;

				case 2:

					sprintf(Problem.Name, "IC_Sawtooth");

					Density_Func_Ptr = &Sawtooth_Density;

					break;

				case 3:

					sprintf(Problem.Name, "IC_SineWave");

					Density_Func_Ptr = &SineWave_Density;

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

			sprintf(Problem.Name, "IC_Magneticum");

			Problem.Boxsize[0] = 1;
			Problem.Boxsize[1] = 1;
			Problem.Boxsize[2] = 0.1;

            Problem.Mpart = 1.0 / Param.Npart * 
				(Problem.Boxsize[0] * Problem.Boxsize[1] * Problem.Boxsize[2]);
	
			Density_Func_Ptr = &Magneticum_Density;

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

void calculateParticleMassByIntegration()
{
    const int N = 100;
    const double dx = Problem.Boxsize[0]/N, dy = Problem.Boxsize[1]/N, dz = Problem.Boxsize[2]/N;
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
}
