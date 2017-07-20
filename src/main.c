#include "globals.h"

struct ParticleData *P;
struct GasParticleData *SphP;

struct Parameters Param = { 0 };
struct ProblemParameters Problem = { 0 };
#ifdef EAT_PNG
struct ImageProperties Image = { 0 };
#endif

float ( *Density_Func_Ptr ) ( const int ) = NULL;
float ( *U_Func_Ptr ) ( const int ) = NULL;
void ( *Magnetic_Field_Func_Ptr ) ( const int, float * ) = NULL;
void ( *Velocity_Func_Ptr ) ( const int, float * ) = NULL;
void ( *PostProcessing_Func_Ptr ) () = NULL;

#pragma omp threadprivate(Omp)
struct OpenMP_infos Omp = { 0 };

/* Set up a box with SPH particles, then WVT relax the particles */

int main ( int argc, char *argv[] )
{
    printf ( "--- This is %s, Version %s ---\n", PROG_NAME, VERSION );

    #pragma omp parallel
    {

        Omp.ThreadID = omp_get_thread_num();
        Omp.NThreads = omp_get_num_threads();
        Omp.Seed[2] = 14041981 * ( Omp.ThreadID + 1 );
        erand48 ( Omp.Seed ); // remove leading 0

        if ( Omp.ThreadID == 0 ) {
            printf ( "Running with %d Threads\n", Omp.NThreads );
#ifdef TWO_DIM
            printf ( "Running in two dimensions\n" );
#else
            printf ( "Running in three dimensions\n" );
#endif
        }

    } // omp parallel

    Assert ( argc == 2, "Usage : ./wvtbox $parameterfile\n" );

#if defined(SPH_CUBIC_SPLINE) && defined(SPH_WC2)
    printf ( "Please decide wether cubic spline or wendland c2 (or standard wc6) shall be used, both may not be definde at the same time\n" );
    return 1;
#endif

    Read_param_file ( argv[1] );

    Setup();

    Make_Positions();

    Make_IDs();

    Regularise_sph_particles();

    Make_Velocities();

    Make_Temperatures();

    Make_Magnetic_Fields();

    Make_PostProcessing();

    Write_output ( 1 ); // verbose

    return EXIT_SUCCESS ;
}

