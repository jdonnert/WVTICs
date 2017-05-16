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

//! @todo  Can be optimized later: can completely forget about saving the positions if we are a bit smartet
double *peanoWalk()
{
    const uint64_t max = peanoCurveLength();
    double *coords = ( double * ) malloc ( 3 * max * sizeof ( double ) );

    uint64_t x, y, z, xmax = 0, ymax = 0, zmax = 0;
    for ( uint64_t peano = 0, i = 0; peano < max; ++peano, i += 3 ) {
        peanoToCoords ( peano, &x, &y, &z );

        coords[i] = x;
        coords[i + 1] = y;
        coords[i + 2] = z;

        //! @todo  Should be able to calculate these directly from the order
        xmax = x > xmax ? x : xmax;
        ymax = y > ymax ? y : ymax;
        zmax = z > zmax ? z : zmax;
    }

    double *cellSides = peanoCellSidelengths ( coords );
    xmax += cellSides[0];
    ymax += cellSides[1];
    zmax += cellSides[2];

    // Translate to center of cell and renormalize box
    for ( uint64_t i = 0; i < 3 * max; i += 3 ) {
        coords[i] = ( coords[i] + 0.5 * cellSides[0] ) * Problem.Boxsize[0] / xmax;
        coords[i + 1] = ( coords[i + 1] + 0.5 * cellSides[1] ) * Problem.Boxsize[1] / ymax;
        coords[i + 2] = ( coords[i + 2] + 0.5 * cellSides[2] ) * Problem.Boxsize[2] / zmax;
    }

    free ( cellSides );

    return coords;
}

//! @todo  Should be able to calculate this directly from the order
double *peanoCellSidelengths ( double *dcoords )
{
    double x1 = dcoords[0], y1 = dcoords[1], z1 = dcoords[2];
    double x2, y2, z2;
    double *volume = ( double * ) malloc ( 3 * sizeof ( double ) );

    const uint64_t max = 1 << ( 3 * PEANO_ORDER );
    for ( uint64_t i = 0; i < 3 * max; i += 3 ) {
        x2 = dcoords[i];
        y2 = dcoords[i + 1];
        z2 = dcoords[i + 2];

        if ( x1 != x2 && y1 != y2 && z1 != z2 ) {
            break;
        }
    }

    volume[0] = x2 - x1;
    volume[1] = y2 - y1;
    volume[2] = z2 - z1;
    return volume;
}

// So far only a test to look at the curve
void peanoWalkToFile()
{
    FILE *pFile;
    pFile = fopen ( "peano.out", "wb" );

    const uint64_t max = peanoCurveLength();
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

#undef PEANO_ORDER
