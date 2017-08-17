#ifndef WVTICS_REDISTRIBUTION_H
#define WVTICS_REDISTRIBUTION_H

#include "globals.h"

void redistributeParticles ( const int movePart );

int findParticleToRedistribute();
int findParticleAsTargetLocation();
void moveParticleInNeighborhoodOf ( const int ipart, const int jpart );

int randomParticle();
bool isOverdense ( const int ipart );
bool isUnderdense ( const int ipart );
float relativeDensityError ( const int ipart );
float getPositionInProximity ( const int jpart, const int i );


#endif //WVTICS_REDISTRIBUTION_H
