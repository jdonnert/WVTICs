#include "globals.h"

void Make_Positions()
{
    printf ( "Sampling positions ...\n" );
    fflush ( stdout );

#ifdef PEANO_SAMPLING
    double *peanoCoords = peanoWalk();
    const uint64_t countCoords = peanoCurveLength();

    const uint64_t cellSize = peanoCellSize();
    const double norm = peanoNormFactor();
    double cellSides[3];
    cellSides[0] = cellSize * Problem.Boxsize[0] * norm;
    cellSides[1] = cellSize * Problem.Boxsize[1] * norm;
    cellSides[2] = cellSize * Problem.Boxsize[2] * norm;

    const double cellVolume = cellSides[0] * cellSides[1] * cellSides[2];

    printf ( " Have %lu peano cells of volume (%g, %g, %g) for %d particles\n", countCoords, cellSides[0], cellSides[1], cellSides[2], Param.Npart );
    Assert ( countCoords > Param.Npart, "Need more peano cells than particles\n" );

    const double probabilityFactor = cellVolume / ( Problem.Mpart );
    double probabilitySum = 0.0;

    int ipart = 0;
    for ( uint64_t i = 0; i < countCoords; ++i ) {
        P[ipart].Pos[0] = peanoCoords[i * 3];
        P[ipart].Pos[1] = peanoCoords[i * 3 + 1];
        P[ipart].Pos[2] = peanoCoords[i * 3 + 2];
        const double probability = probabilityFactor * Density_Func_Ptr ( ipart );
        probabilitySum += probability;

        //Accept particle
        if ( ipart < Param.Npart && probability > erand48 ( Omp.Seed ) ) {
            //randomize position inside peano cell
            P[ipart].Pos[0] += ( erand48 ( Omp.Seed ) - 0.5 ) * cellSides[0];
            P[ipart].Pos[1] += ( erand48 ( Omp.Seed ) - 0.5 ) * cellSides[1];
            P[ipart].Pos[2] += ( erand48 ( Omp.Seed ) - 0.5 ) * cellSides[2];

            ++ipart;

            if ( ipart == Param.Npart ) {
                printf ( " Aborting at %lu of %lu peano nodes (%g%%)\n", i, countCoords, i * 100. / countCoords );
                //Let it run further to get correct probabilitySum values
            }
        }
    }
    if ( ipart != Param.Npart ) {
        Problem.Mpart = Problem.Mpart * ipart / Param.Npart;
        Param.Npart = ipart; //If got less particles we can just ignore all memory beyond this point in the structs
        printf ( " Resetting particle number to %d and particle mass to %g\n", Param.Npart, Problem.Mpart );
    }

    // Normalization: Sum_cells p = 1 * Npart
    printf ( " Sum of probabilities %g (%g%%)\n", probabilitySum, probabilitySum * 100.0 / Param.Npart );

    free ( peanoCoords );
#else
    #pragma omp parallel for
    for ( int ipart = 0; ipart < Param.Npart; ipart++ ) {
#ifdef REJECTION_SAMPLING
        double rho = 0.0, rho_r = 0.0;

        while ( rho >= rho_r ) {

            P[ipart].Pos[0] = erand48 ( Omp.Seed ) * Problem.Boxsize[0];
            P[ipart].Pos[1] = erand48 ( Omp.Seed ) * Problem.Boxsize[1];
            P[ipart].Pos[2] = erand48 ( Omp.Seed ) * Problem.Boxsize[2];

            rho = Problem.Rho_Max * erand48 ( Omp.Seed );
            rho_r = Density_Func_Ptr ( ipart );
        }
#else
        P[ipart].Pos[0] = erand48 ( Omp.Seed ) * Problem.Boxsize[0];
        P[ipart].Pos[1] = erand48 ( Omp.Seed ) * Problem.Boxsize[1];
        P[ipart].Pos[2] = erand48 ( Omp.Seed ) * Problem.Boxsize[2];
#endif //REJECTION_SAMPLING

        P[ipart].Type = 0;
    }
#endif //PEANO_SAMPLING

    printf ( "done\n" );

    return;
}

void Make_Velocities()
{

    printf ( "Velocities ..." );
    fflush ( stdout );

    #pragma omp parallel for
    for ( int ipart = 0; ipart < Param.Npart; ipart++ ) {

        ( *Velocity_Func_Ptr ) ( ipart, P[ipart].Vel );
    }

    printf ( " done\n" );

    return ;
}

void Make_Temperatures()
{

    printf ( "Internal Energy ..." );
    fflush ( stdout );

    #pragma omp parallel for
    for ( int ipart = 0; ipart < Param.Npart; ipart++ ) {

        SphP[ipart].U = ( *U_Func_Ptr ) ( ipart );
    }

    printf ( " done\n" );

    return ;
}

void Make_Magnetic_Fields()
{

    printf ( "Magnetic Field ..." );
    fflush ( stdout );

    #pragma omp parallel for
    for ( int ipart = 0; ipart < Param.Npart; ipart++ ) {

        ( *Magnetic_Field_Func_Ptr ) ( ipart, SphP[ipart].Bfld );

    }

    printf ( " done\n" );

    return ;
}

void Make_PostProcessing()
{

    printf ( "Post Processing ..." );
    fflush ( stdout );

    ( *PostProcessing_Func_Ptr ) ();

    printf ( " done\n" );

    return ;
}

