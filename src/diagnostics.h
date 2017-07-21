#ifndef WVTICS_DIAGNOSTICS_H
#define WVTICS_DIAGNOSTICS_H

struct Quadruplet {
    double min, max, mean, sigma;
};

struct Quadruplet calculateStatsOn ( float *values[3], int n );

void writeIterationDiagnostics ( const int iteration,
                                 const struct Quadruplet *error,
                                 const double diffError,
                                 const double moveMps[4],
                                 const struct Quadruplet *delta );

#endif //WVTICS_DIAGNOSTICS_H
