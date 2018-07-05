#ifndef WVTICS_PEANOWALK_H
#define WVTICS_PEANOWALK_H

#include "external/morton_utils.h"

#define PEANO_ORDER 21
typedef int64_t peanoKey;

void peanoToCoords ( uint64_t peano, uint64_t *x, uint64_t *y, uint64_t *z );
uint64_t coordsToPeano ( uint64_t x, uint64_t y, uint64_t z );

uint64_t peanoCurveLength();
uint64_t peanoMaxCoord();
uint64_t peanoCellSize();
double peanoNormFactor();

void assignPeanoCoordinates ( float *coords, uint64_t peano );
void translateAndRenormalizePeanoCoords ( float *coords, const double halfCellSize, const double norm );
void translateAndRenormalizeCoordsForPeano ( float *coords, const double halfCellSize, const double norm );

void Print ( const __uint128_t val );
void Printr ( const __uint128_t val );

/*unsigned char reverseBitsInByte(unsigned char b);
uint32_t reverseBitsInInt32(uint32_t i);
uint64_t reverseBitsInInt64(uint64_t i);*/

uint64_t reverseBitsInPeanoKey ( uint64_t key );


void Sort_Particles_By_Peano_Key();
peanoKey Peano_Key ( const double, const double, const double );
peanoKey Reversed_Peano_Key ( const double, const double, const double );
void test_peanokey();

#endif //WVTICS_PEANOWALK_H
