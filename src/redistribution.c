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

//! @todo count rejected particles and set a max for that; maybe we can remove the decay of the redistribution then again
void redistributeParticles ( const int movePart )
{
    printf ( "Redistributing %d particles (=%g%%)\n", movePart, movePart * 100. / Param.Npart );

    resetRedistributionFlags();

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
        while ( P[ipart].Redistributed || !acceptParticleForMovement ( ipart ) ) {
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

    while ( !acceptParticleAsTarget ( ipart ) ) {
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

double randomNumber()
{
    return erand48 ( Omp.Seed );
}

int randomParticle()
{
    return randomNumber() * Param.Npart;
}

bool acceptParticleForMovement ( const int ipart )
{
    const double r = randomNumber();
    return r < erf ( relativeDensityErrorWithSign ( ipart ) );
}

bool acceptParticleAsTarget ( const int ipart )
{
    const double r = randomNumber();
    return r < ( -1.0 ) * relativeDensityErrorWithSign ( ipart );
}

//This ranges from -1 (Rho = 0) to infinity (Rho = infinity)
float relativeDensityErrorWithSign ( const int ipart )
{
    const float rhoModel = ( *Density_Func_Ptr ) ( ipart );
    return ( SphP[ipart].Rho - rhoModel ) / rhoModel;
}

float relativeDensityError ( const int ipart )
{
    const float rho = ( *Density_Func_Ptr ) ( ipart );
    return fabs ( relativeDensityErrorWithSign ( ipart ) );
}

float getPositionInProximity ( const int ipart, const int i )
{
    float ret = -1.0;

    while ( ret < 0.0 || ret >= Problem.Boxsize[i] ) {
        ret = P[ipart].Pos[i] + ( 2.0 * erand48 ( Omp.Seed ) - 1.0 ) * SphP[ipart].Hsml * 0.3;
        if ( Problem.Periodic ) {
            if ( ret >= Problem.Boxsize[i] ) {
                ret -= Problem.Boxsize[i];
            } else if ( ret < 0.0 ) {
                ret += Problem.Boxsize[i];
            }
        }
    }

    return ret;
}

void resetRedistributionFlags()
{
    #pragma omp parallel for
    for ( int ipart = 0; ipart < Param.Npart; ++ipart ) {
        P[ipart].Redistributed = false;
    }
}
