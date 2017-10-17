#ifndef WVTICS_REDISTRIBUTION_H
#define WVTICS_REDISTRIBUTION_H

#include "globals.h"

void redistributeParticles ( const int movePart, const int maxProbes );

int findParticleToRedistribute ( int *probes );
int findParticleAsTargetLocation();
void moveParticleInNeighborhoodOf ( const int ipart, const int jpart );

double randomNumber();
int randomParticle();

bool acceptParticleForMovement ( const int ipart );
bool acceptParticleAsTarget ( const int ipart );
float relativeDensityErrorWithSign ( const int ipart );
float relativeDensityError ( const int ipart );
float getPositionInProximity ( const int ipart, const int i );

void resetRedistributionFlags();


#endif //WVTICS_REDISTRIBUTION_H
