#ifndef WVTICS_REDISTRIBUTION_H
#define WVTICS_REDISTRIBUTION_H

#include "globals.h"

void redistributeParticles();

int findParticleToRedistribute();
int findParticleAsTargetLocation();
void moveParticleInNeighborhoodOf ( const int ipart, const int jpart );

int randomParticle();
bool isOverdense ( const int ipart );
float relativeDensityError ( const int ipart );


#endif //WVTICS_REDISTRIBUTION_H
