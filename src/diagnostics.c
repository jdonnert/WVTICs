#include "diagnostics.h"
#include "globals.h"
#include <float.h>
#include <stdio.h>

struct Quadruplet calculateStatsOn ( float *values[3], int n )
{
    struct Quadruplet q;
    double min = DBL_MAX, max = 0., mean = 0., sigma = 0.;

    #pragma omp parallel for reduction(+:mean,sigma) reduction(max:max) reduction(min:min)
    for ( int  i = 0; i < n; i++ ) {
        const float v2 = p2 ( values[0][i] ) + p2 ( values[1][i] ) + p2 ( values[2][i] );
        const float v = sqrt ( v2 );

        min = fmin ( v, min );
        max = fmax ( v, max );
        mean += v;
        sigma += v2;
    }
    mean /= n;
    sigma = sqrt ( sigma - p2 ( mean ) );

    q.min = min;
    q.max = max;
    q.mean = mean;
    q.sigma = sigma;

    return q;
}


void initIterationDiagnostics()
{
    FILE *fp;

    if ( ! ( fp = fopen ( "diagnostics.log", "w" ) ) ) {
        fprintf ( stderr, "Can't open file diagnostics.log\n" );
    }

    fprintf ( fp, "Iter\tError min\tError max\tError mean\tError sigma\tError diff\tMove dmps\tdmps/10\t\
            dmps/100\tdmps/1000\tDelta min\tDelta max\tDelta mean\tDelta sigma\n" );

    fclose ( fp );
}


void writeIterationDiagnostics ( const int iteration,
                                 const struct Quadruplet *error,
                                 const double diffError,
                                 const double moveMps[4],
                                 const struct Quadruplet *delta )
{
    FILE *fp;

    if ( ! ( fp = fopen ( "diagnostics.log", "a" ) ) ) {
        fprintf ( stderr, "Can't open file diagnostics.log\n" );
    }

    fprintf ( fp, "%03d\t%+7.5e\t%+7.5e\t%+7.5e\t%+7.5e\t%+7.5e\t%+7.5e\t%+7.5e\t%+7.5e\t%+7.5e\t%+7.5e\t%+7.5e\t%+7.5e\t%+7.5e\n",
              iteration,
              error->min,
              error->max,
              error->mean,
              error->sigma,
              diffError,
              moveMps[0],
              moveMps[1],
              moveMps[2],
              moveMps[3],
              delta->min,
              delta->max,
              delta->mean,
              delta->sigma );

    fclose ( fp );
}
