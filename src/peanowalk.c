#include <printf.h>
#include <stdio.h>
#include <stdlib.h>
#include "peanowalk.h"
#include "globals.h"

#define PEANO_ORDER 9

void peanoToCoords ( uint64_t peano, uint64_t *x, uint64_t *y, uint64_t *z )
{
    //Peano -> Morton
    const int bits = 21;
    uint64_t morton = peano2morton3D_64 ( peano, bits );

    //Morton -> Coords
    morton2coord3D_64 ( morton, x, y, z );
}

uint64_t peanoCurveLength()
{
    return 1 << ( 3 * PEANO_ORDER );
}

uint64_t peanoMaxCoord()
{
    return ( 1 << PEANO_ORDER ) - 1;
}

double peanoNormFactor()
{
    return 1.0 / ( peanoMaxCoord() + peanoCellSize() );
}

uint64_t peanoCellSize()
{
    return 1;
}

//! @todo  Can be optimized later: can completely forget about saving the positions if we are a bit smartet
double *peanoWalk()
{
    const uint64_t max = peanoCurveLength();
    double *coords = ( double * ) malloc ( 3 * max * sizeof ( double ) );

    uint64_t x, y, z;
    for ( uint64_t peano = 0, i = 0; peano < max; ++peano, i += 3 ) {
        peanoToCoords ( peano, &x, &y, &z );

        coords[i] = x;
        coords[i + 1] = y;
        coords[i + 2] = z;
    }

    const uint64_t cellSize = peanoCellSize();
    const double halfCellSize = 0.5 * cellSize;

    const double norm = peanoNormFactor();

    // Translate to center of cell and renormalize box
    for ( uint64_t i = 0; i < 3 * max; i += 3 ) {
        coords[i] = ( coords[i] + halfCellSize ) * Problem.Boxsize[0] * norm;
        coords[i + 1] = ( coords[i + 1] + halfCellSize ) * Problem.Boxsize[1] * norm;
        coords[i + 2] = ( coords[i + 2] + halfCellSize ) * Problem.Boxsize[2] * norm;
    }

    return coords;
}

#undef PEANO_ORDER
