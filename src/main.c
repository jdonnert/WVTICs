#include "globals.h"

struct ParticleData *P;
struct GasParticleData *SphP;
struct Units Unit;

struct Parameters Param = { 0 };
struct ProblemParameters Problem = { 0 };

float (*Density_Func_Ptr) (const int) = NULL;
float (*U_Func_Ptr) (const int) = NULL;
void (*Magnetic_Field_Func_Ptr) (const int, float *) = NULL;
void (*Velocity_Func_Ptr) (const int, float *) = NULL;

#pragma omp threadprivate(Omp)
struct OpenMP_infos Omp = { 0 };

/* Set up a box with SPH particles, then WVT relax the particles */

int main(int argc, char *argv[])
{
    printf("--- This is %s, Version %s ---\n", PROG_NAME, VERSION);

#pragma omp parallel
    {

        Omp.ThreadID = omp_get_thread_num();
        Omp.NThreads = omp_get_num_threads();
        Omp.Seed[2] = 14041981 * (Omp.ThreadID + 1);
        erand48(Omp.Seed); // remove leading 0

        if (Omp.ThreadID == 0)
            printf("Running with %d Threads\n", Omp.NThreads);

    } // omp parallel

    Assert(argc == 2, "Usage : ./wvtbox $parameterfile\n");

    Read_param_file(argv[1]);

    Set_units();

    Setup();

    Make_positions();

    Make_IDs();
    
    Regularise_sph_particles();

 	Make_Velocities();

	Make_Temperatures();
	
	Make_Magnetic_Field();

    Write_output();

    return EXIT_SUCCESS ;
}

