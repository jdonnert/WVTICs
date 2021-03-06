#ifndef WVTICS_PEANOWALK_H
#define WVTICS_PEANOWALK_H

#include "external/morton_utils.h"

void peanoToCoords ( uint64_t peano, uint64_t *x, uint64_t *y, uint64_t *z );

uint64_t peanoCurveLength();
uint64_t peanoMaxCoord();
uint64_t peanoCellSize();
double peanoNormFactor();

void assignPeanoCoordinates ( float *coords, uint64_t peano );
void translateAndRenormalizePeanoCoords ( float *coords, const double halfCellSize, const double norm );

#endif //WVTICS_PEANOWALK_H
