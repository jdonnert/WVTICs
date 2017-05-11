#ifndef WVTICS_PEANOWALK_H
#define WVTICS_PEANOWALK_H

#include "external/morton_utils.h"

void peanoToCoords ( uint64_t peano, uint64_t *x, uint64_t *y, uint64_t *z );

void peanoWalk();

#endif //WVTICS_PEANOWALK_H
