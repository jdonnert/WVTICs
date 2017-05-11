#include <printf.h>
#include <stdio.h>
#include "peanowalk.h"

void peanoToCoords ( uint64_t peano, uint64_t *x, uint64_t *y, uint64_t *z )
{
    //Peano -> Morton
    const int bits = 21;
    uint64_t morton = peano2morton3D_64 ( peano, bits );

    //Morton -> Coords
    morton2coord3D_64 ( morton, x, y, z );
}

// So far only a test to look at the curve
void peanoWalk()
{
    FILE *pFile;
    pFile = fopen ( "peano.out", "wb" );

    const int n = 4;
    const uint64_t max = 1 << ( 3 * n );
    printf ( "This is going to be %g GB (%lu points at curve)\n", max * 3.0 * 8 / 1024 / 1024 / 1024, max );

    uint64_t x, y, z;
    for ( uint64_t peano = 0; peano < max; ++peano ) {
        peanoToCoords ( peano, &x, &y, &z );
        fwrite ( &x, sizeof ( uint64_t ), 1, pFile );
        fwrite ( &y, sizeof ( uint64_t ), 1, pFile );
        fwrite ( &z, sizeof ( uint64_t ), 1, pFile );
    }

    fclose ( pFile );
}
