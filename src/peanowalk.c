#include <printf.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_heapsort.h>
#include "peanowalk.h"
#include "globals.h"

#define MORTON_BITS 21

void peanoToCoords ( uint64_t peano, uint64_t *x, uint64_t *y, uint64_t *z )
{
    //Peano -> Morton
    uint64_t morton = peano2morton3D_64 ( peano, MORTON_BITS );

    //Morton -> Coords
    morton2coord3D_64 ( morton, x, y, z );
}

uint64_t coordsToPeano ( uint64_t x, uint64_t y, uint64_t z )
{
    //Coords -> Morton
    uint64_t morton = coord2morton3D_64 ( x, y, z );

    //Mortron -> Peano
    uint64_t peano = morton2peano3D_64 ( morton, MORTON_BITS );
    return peano;
}

uint64_t peanoCurveLength()
{
    return 1LL << ( 3LL * PEANO_ORDER );
}

uint64_t peanoMaxCoord()
{
    return ( 1LL << PEANO_ORDER ) - 1LL;
}

double peanoNormFactor()
{
    return 1.0 / ( peanoMaxCoord() + peanoCellSize() );
}

uint64_t peanoCellSize()
{
    return 1;
}

void assignPeanoCoordinates ( float *coords, uint64_t peano )
{
    uint64_t x, y, z;
    peanoToCoords ( peano, &x, &y, &z );

    coords[0] = x;
    coords[1] = y;
    coords[2] = z;
}

void translateAndRenormalizePeanoCoords ( float *coords, const double halfCellSize, const double norm )
{
    coords[0] = ( coords[0] + halfCellSize ) * Problem.Boxsize[0] * norm;
    coords[1] = ( coords[1] + halfCellSize ) * Problem.Boxsize[1] * norm;
    coords[2] = ( coords[2] + halfCellSize ) * Problem.Boxsize[2] * norm;
}

void translateAndRenormalizeCoordsForPeano ( float *coords, const double halfCellSize, const double norm )
{
    coords[0] = coords[0] / ( Problem.Boxsize[0] * norm ); // - halfCellSize;
    coords[1] = coords[1] / ( Problem.Boxsize[1] * norm ); // - halfCellSize;
    coords[2] = coords[2] / ( Problem.Boxsize[2] * norm ); // - halfCellSize;
}

#ifdef NEW_PEANO

void Print ( const __uint128_t val )
{
    for ( int i = 63; i >= 0; i-- ) {

        printf ( "%llu", ( long long ) ( ( val & ( ( __uint128_t ) 1 << i ) ) >> i ) );

        if ( i % 3 - 1 == 0 && i != 0 ) {
            printf ( "." );
        }
    }
    printf ( "\n" );
    fflush ( stdout );

    return ;
}


void Printr ( const __uint128_t val )
{
    for ( int i = 63; i >= 0; i-- ) {

        printf ( "%llu", ( long long ) ( ( val & ( ( __uint128_t ) 1 << i ) ) >> i ) );

        if ( i % 3 == 0 && i != 0 ) {
            printf ( "." );
        }
    }
    printf ( "\n" );
    fflush ( stdout );

    return ;
}

uint64_t reverseBitsInPeanoKey ( uint64_t key )
{
    //printf("Reverse...\n");
    uint64_t reversed = key & ( 0x1LL << ( 3LL * PEANO_ORDER ) );
    //printf("Starting: ");
    //Printr(reversed);

    for ( int whichTriplet = 0; whichTriplet < PEANO_ORDER; ++ whichTriplet ) { // breaks with < 11
        const uint64_t skip1 = 3 * whichTriplet;
        const uint64_t skip2 = 3 * ( PEANO_ORDER - whichTriplet - 1 ) + 1;

        const uint64_t inTriplet = ( key & ( 0x7LL << skip1 ) ) >> skip1;
        //printf("Selector: ");
        //Print(inTriplet);
        reversed |= inTriplet << skip2;
        //printf("Result:   ");
        //Printr(reversed);
    }

    return reversed;
}



int compare_peanoKeys ( const void *a, const void *b )
{
    const peanoKey *x = ( const peanoKey * ) a;
    const peanoKey *y = ( const peanoKey * ) b;

    return ( int ) ( *x > *y ) - ( *x < *y );
}

static peanoKey *Keys = NULL;
static size_t *Idx = NULL;

static void reorder_particles()
{
    for ( int i = 0; i < Param.Npart; i++ ) {

        if ( Idx[i] == i ) {
            continue;
        }

        int dest = i;

        struct ParticleData Ptmp = P[i];
        struct GasParticleData Sphtmp = SphP[i];

        int src = Idx[i];

        for ( ;; ) {

            P[dest] = P[src];
            SphP[dest] = SphP[src];

            Idx[dest] = dest;

            dest = src;

            src = Idx[dest];

            if ( src == i ) {
                break;
            }
        }

        P[dest] = Ptmp;
        SphP[dest] = Sphtmp;

        Idx[dest] = dest;

    } // for i

    return ;
}

void Sort_Particles_By_Peano_Key()
{
    if ( Keys == NULL ) {
        Keys = malloc ( Param.Npart * sizeof ( *Keys ) );
    } else {
        memset ( Keys, 0, Param.Npart * sizeof ( *Keys ) );
    }

    if ( Idx == NULL ) {
        Idx = malloc ( Param.Npart * sizeof ( *Idx ) );
    } else {
        memset ( Idx, 0, Param.Npart * sizeof ( *Idx ) );
    }

    #pragma omp parallel for
    for ( int ipart = 0; ipart < Param.Npart; ipart++ ) {

        double px = P[ipart].Pos[0] / Problem.Boxsize[0];
        double py = P[ipart].Pos[1] / Problem.Boxsize[1];
        double pz = P[ipart].Pos[2] / Problem.Boxsize[2];

        P[ipart].Key = Keys[ipart] = Peano_Key ( px, py, pz );
    }

    gsl_heapsort_index ( Idx, Keys, Param.Npart, sizeof ( *Keys ), &compare_peanoKeys );

    reorder_particles();

    return ;
}


peanoKey Peano_Key ( const double x, const double y, const double z )
{
    // @todo pull these out for performance gain
    const uint64_t cellSize = peanoCellSize();
    const double halfCellSize = 0.5 * cellSize;
    const double norm = peanoNormFactor();

    float a[3] = {x, y, z};

    translateAndRenormalizeCoordsForPeano ( a, halfCellSize, norm );
    int64_t key = coordsToPeano ( a[0], a[1], a[2] );
    key = key << 0x1; // useless bit from start to end

    return key;
}

peanoKey Reversed_Peano_Key ( const double x, const double y, const double z )
{
    // @todo pull these out for performance gain
    const uint64_t cellSize = peanoCellSize();
    const double halfCellSize = 0.5 * cellSize;
    const double norm = peanoNormFactor();

    float a[3] = {x, y, z};

    translateAndRenormalizeCoordsForPeano ( a, halfCellSize, norm );
    int64_t key = coordsToPeano ( a[0], a[1], a[2] );
    key = reverseBitsInPeanoKey ( key );
    key = ( key >> 0x1 ) & ~ ( 1LL << 63 ); // useless bit from end to start
    key = ( key << 3 ); // Patch: add root node on the right

    return key;
}

void test_peanokey()
{
    Problem.Boxsize[0] = 1.0;
    Problem.Boxsize[1] = 1.0;
    Problem.Boxsize[2] = 1.0;

    float a[3] = { 0 };
    for ( int order = 0; order < 2; ++order ) {
        printf ( "\nTesting order %d\n\n", order );
        float delta = 1 / pow ( 2.0, order );
        int n = roundf ( 1 / delta );

        const uint64_t cellSize = peanoCellSize();
        const double halfCellSize = 0.5 * cellSize;
        const double norm = peanoNormFactor();
        printf ( "Half cell size: %g; norm: %g\n", halfCellSize, norm );

        for ( int i = 0; i < n; i++ ) {
            for ( int j = 0; j < n; j++ ) {
                for ( int k = 0; k < n; k++ ) {

                    a[0] = ( i + 0.5 ) * delta / Problem.Boxsize[0];
                    a[1] = ( j + 0.5 ) * delta / Problem.Boxsize[1];
                    a[2] = ( k + 0.5 ) * delta / Problem.Boxsize[2];

                    peanoKey stdkey = Peano_Key ( a[0], a[1], a[2] );

                    printf ( "%g %g %g %llu  \n", a[0], a[1], a[2],
                             ( long long ) stdkey );

                    Print ( stdkey );

                    peanoKey reversed = Reversed_Peano_Key ( a[0], a[1], a[2] );

                    Printr ( reversed );

                    printf ( "\n" );
                }
            }
        }
    }

    return ;
}

#endif

#undef MORTON_BITS
