#include "redistribution.h"

/*
 * This works similar to the Metropolis algorithm:
 * We find a random untouched particle with density higher than the model (= high energy state)
 * With a small probability we also accept lower densities in order to break out of local minima (to be implemented)
 * We find a random underdense region (= low energy state)
 * We move the particle in the proximity of that region
 *
 * We do that for a small percentage of particles and assume that this does not change the density to much,
 * so that we can omit recalculation of density during the process.
 */

void redistributeParticles ( const int movePart )
{
    printf ( "Redistributing %d particles (=%g%%)\n", movePart, movePart * 100. / Param.Npart );

    #pragma omp parallel for
    for ( int i = 0; i < movePart; ++i ) {
        const int ipart = findParticleToRedistribute();
        const int jpart = findParticleAsTargetLocation();
        moveParticleInNeighborhoodOf ( ipart, jpart );
    }
}

int findParticleToRedistribute()
{
    int ipart = randomParticle();

    //! @todo eventually we want to even accept underdense particles with a small probability
    bool run = true;
    while ( run ) {
        while ( P[ipart].Redistributed || isUnderdense ( ipart ) ) {
            ipart = randomParticle();
        }
        #pragma omp critical
        if ( !P[ipart].Redistributed ) {
            P[ipart].Redistributed = true;
            run = false;
        }
    }

    return ipart;
}

int findParticleAsTargetLocation()
{
    int ipart = randomParticle();

    while ( isOverdense ( ipart ) ) {
        ipart = randomParticle();
    }

    return ipart;
}

void moveParticleInNeighborhoodOf ( const int ipart, const int jpart )
{
    for ( int i = 0; i < 3; ++i ) {
        P[ipart].Pos[i] = getPositionInProximity ( jpart, i );
    }
}

int randomParticle()
{
    return erand48 ( Omp.Seed ) * Param.Npart;
}

bool isOverdense ( const int ipart )
{
    return SphP[ipart].Rho > ( *Density_Func_Ptr ) ( ipart );
}

bool isUnderdense ( const int ipart )
{
    return SphP[ipart].Rho < ( *Density_Func_Ptr ) ( ipart );
}

float relativeDensityError ( const int ipart )
{
    const float rho = ( *Density_Func_Ptr ) ( ipart );
    return fabs ( SphP[ipart].Rho - rho ) / rho;
}

float getPositionInProximity ( const int jpart, const int i )
{
    float ret = -1.0;

    while ( ret < 0.0 || ret >= Problem.Boxsize[i] ) {
        ret = P[jpart].Pos[i] + erand48 ( Omp.Seed ) * SphP[jpart].Hsml * 0.3;
    }

    return ret;
}
