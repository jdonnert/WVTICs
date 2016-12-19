#include "globals.h"

float zero_function(const int ipart);
void zero_function_vec(const int ipart, float out[3]);

void setup_problem(const int Flag, const int Subflag);

void Setup()
{
	setup_problem(Param.Problem_Flag, Param.Problem_Subflag);

	printf( "Problem %d.%d : %s \n"
			"   Npart: %d \n"
			"   Mpart: %g \n"
			"   Boxsize: %g x %g x %g \n"
            "   Periodic: %d \n\n"
			,Param.Problem_Flag, Param.Problem_Subflag, Problem.Name,
			Param.Npart, Problem.Mpart, Problem.Boxsize[0], 
			Problem.Boxsize[1], Problem.Boxsize[2], Problem.Periodic );

    size_t nBytes = Param.Npart * sizeof(*P); // allocate particles

    P = Malloc(nBytes);
    
	memset(P, 0, nBytes);

    nBytes = Param.Npart * sizeof(*SphP);
    
	SphP = Malloc(nBytes);
    
	memset(SphP, 0, nBytes);

    return;
}

void setup_problem(const int Flag, const int Subflag)
{
	Density_Func_Ptr = &zero_function;
	U_Func_Ptr = &zero_function;
	Velocity_Func_Ptr = &zero_function_vec;
	Magnetic_Field_Func_Ptr = &zero_function_vec;

    Problem.Periodic = true;
	Problem.Rho_Max = 1.0;

	switch (Flag) {
		case 0:

			switch (Subflag) {
				case 0:

					sprintf(Problem.Name, "IC_Constant_Density");

					Problem.Boxsize[0] = Problem.Boxsize[1] = Problem.Boxsize[2] = 1;
					Problem.Mpart = 1.0 / Param.Npart;

					Density_Func_Ptr = &Constant_Density;

					break;

				case 1:

					sprintf(Problem.Name, "IC_TopHat");

					Problem.Boxsize[0] = Problem.Boxsize[1] = Problem.Boxsize[2] = 1;
					Problem.Mpart = 1.0 / Param.Npart;

					Density_Func_Ptr = &TopHat_Density;

					break;

				case 2:

					sprintf(Problem.Name, "IC_Sawtooth");

					Problem.Boxsize[0] = Problem.Boxsize[1] = Problem.Boxsize[2] = 1;
					Problem.Mpart = 1.0 / Param.Npart;

					Density_Func_Ptr = &Sawtooth_Density;

					break;

				case 3:

					sprintf(Problem.Name, "IC_SineWave");

					Problem.Boxsize[0] = Problem.Boxsize[1] = Problem.Boxsize[2] = 1;
					Problem.Mpart = 1.0 / Param.Npart;

					Density_Func_Ptr = &SineWave_Density;

					break;

				default:

					Assert(false, "Effect %d.%d not implemented", Flag, Subflag);

					break;
			}
			break;
        case 1:

            Problem.Periodic = false;

            switch (Subflag) {
                case 0:

                    sprintf(Problem.Name, "IC_GradientDensity");

                    Problem.Boxsize[0] = Problem.Boxsize[1] = Problem.Boxsize[2] = 1;
                    Problem.Mpart = 1.0 / Param.Npart;

                    Density_Func_Ptr = &Gradient_Density;

                    break;

                default:

                    Assert(false, "Effect %d.%d not implemented", Flag, Subflag);

                    break;
            }
            break;

		case 2:

			sprintf(Problem.Name, "IC_Magneticum");

			Problem.Periodic = false;
			Problem.Boxsize[0] = Problem.Boxsize[1] = Problem.Boxsize[2] = 2;
    		Problem.Mpart = 1.0 / Param.Npart;
	
			Density_Func_Ptr = &Magneticum_Density;

			break;

		default:

			Assert(false, "Effect %d.%d not implemented", Flag, Subflag);

			break;
	}

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
