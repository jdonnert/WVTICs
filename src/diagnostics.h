#ifndef WVTICS_DIAGNOSTICS_H
#define WVTICS_DIAGNOSTICS_H

struct Diagnostics {
    FILE *file;
};

struct Quadruplet {
    double min, max, mean, sigma;
};

Quadruplet calculateStatsOn ( const float *values[3], int n );

void writeIterationDiagnostics ( const int iteration,
                                 const Quadruplet *error,
                                 const double diffError,
                                 const double moveMps[4],
                                 const Quadruplet *delta );

#endif //WVTICS_DIAGNOSTICS_H
