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

#undef PEANO_ORDER
