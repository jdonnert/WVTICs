#include "globals.h"
#include "tree.h"
#include "kernel.h"

#ifdef OUTPUT_DIAGNOSTICS
#include "diagnostics.h"
#endif

#define WVTNNGB DESNNGB

void writeStepFile ( int it );

/* Settle SPH particle with weighted Voronoi tesselations (Diehl+ 2012).
 * Here hsml is not the SPH smoothing length, but is related to a local
 * metric defined ultimately by the density model.
 * Relaxation is done in units of the boxsize, hence the box volume is 1
 * Return code true means that rerun could be usefull */

void Regularise_sph_particles()
{
    const int nPart = Param.Npart;

    const double boxsize[3] = { Problem.Boxsize[0], Problem.Boxsize[1],
                                Problem.Boxsize[2]
                              };
    const double boxhalf[3] = { boxsize[0] / 2, boxsize[1] / 2, boxsize[2] / 2, };

    const double median_boxsize = fmax ( boxsize[1], boxsize[2] ); // boxsize[0] is largest

    printf ( "Starting iterative SPH regularisation \n"
             "   Maxiter=%d, MpsFraction=%g StepReduction=%g LimitMps=(%g,%g,%g,%g)\n\n",
             Param.Maxiter, Param.MpsFraction, Param.StepReduction, Param.LimitMps[0], Param.LimitMps[1], Param.LimitMps[2], Param.LimitMps[3] );
    fflush ( stdout );

    float *hsml = NULL;
    size_t nBytes = nPart * sizeof ( *hsml );
    hsml = Malloc ( nBytes );

    float *delta[3] = { NULL };

    nBytes = nPart * sizeof ( **delta );

    delta[0] = Malloc ( nBytes );
    delta[1] = Malloc ( nBytes );
    delta[2] = Malloc ( nBytes );

    int it = 0;

#ifdef TWO_DIM
    double npart_1D = pow ( nPart, 1.0 / 2.0 );
#else
    double npart_1D = pow ( nPart, 1.0 / 3.0 );
#endif

    double step = 1.0 / ( npart_1D * Param.MpsFraction );

    double errLast = DBL_MAX;
    double errDiff = DBL_MAX;

    const double volume = Problem.Boxsize[0] * Problem.Boxsize[1] * Problem.Boxsize[2];
    const double rho_mean = nPart * Problem.Mpart / volume;

    double last_cnt = DBL_MAX;

    for ( ;; ) {

        Find_sph_quantities();

        if ( it++ > Param.Maxiter ) {
            break;
        }

#ifdef SAVE_WVT_STEPS
        writeStepFile ( it );
#endif

        double errMin = DBL_MAX, errMax = 0, errMean = 0, errSigma = 0.;

        #pragma omp parallel for reduction(+:errMean,errSigma) reduction(max:errMax) reduction(min:errMin)
        for ( int  ipart = 0; ipart < nPart; ipart++ ) { // get error

            float rho = ( *Density_Func_Ptr ) ( ipart );

            float err = fabs ( SphP[ipart].Rho - rho ) / rho;

            errMin = fmin ( err, errMin );
            errMax = fmax ( err, errMax );

            errMean += err;
            errSigma += err * err;
        }

        errMean /= nPart;
        errSigma = sqrt ( errSigma - p2 ( errMean ) );

        errDiff = ( errLast - errMean ) / errMean;

        printf ( "   #%02d: Err min=%3g max=%3g mean=%03g sigma=%03g diff=%03g step=%g\n", it, errMin, errMax, errMean, errSigma, errDiff, step );

        errLast = errMean;

        double vSphSum = 0; // total volume defined by hsml
        double max_hsml = 0;

        #pragma omp parallel for shared(hsml) reduction(+:vSphSum,max_hsml)
        for ( int ipart = 0; ipart < nPart; ipart++ ) { // find hsml

            float rho = ( *Density_Func_Ptr ) ( ipart );

            SphP[ipart].Rho_Model = rho;

#ifdef TWO_DIM
            hsml[ipart] = pow ( WVTNNGB * Problem.Mpart / rho / pi, 1. / 2. );
            vSphSum += p2 ( hsml[ipart] );
#else
            hsml[ipart] = pow ( WVTNNGB * Problem.Mpart / rho / fourpithird, 1. / 3. );
            vSphSum += p3 ( hsml[ipart] );
#endif

            max_hsml = max ( max_hsml, hsml[ipart] );
        }

#ifdef TWO_DIM
        float norm_hsml = pow ( WVTNNGB / vSphSum / pi , 1.0 / 2.0 ) * median_boxsize;
#else
        float norm_hsml = pow ( WVTNNGB / vSphSum / fourpithird , 1.0 / 3.0 ) * median_boxsize;
#endif

        #pragma omp parallel for
        for ( int ipart = 0; ipart < nPart; ipart++ ) {
            hsml[ipart] *= norm_hsml;
        }

        #pragma omp parallel for shared(delta, hsml, P) schedule(dynamic, nPart/Omp.NThreads/256)
        for ( int ipart = 0; ipart < nPart; ipart++ ) {

            delta[0][ipart] = delta[1][ipart] = delta[2][ipart] = 0;

            int ngblist[NGBMAX] = { 0 };
            int ngbcnt = Find_ngb ( ipart, hsml[ipart], ngblist );

#ifdef OUTPUT_DIAGNOSTICS
#endif
            for ( int i = 0; i < ngbcnt; i++ ) { // neighbour loop

                int jpart = ngblist[i];

                if ( ipart == jpart ) {
                    continue;
                }

                float dx = P[ipart].Pos[0] - P[jpart].Pos[0];
                float dy = P[ipart].Pos[1] - P[jpart].Pos[1];
                float dz = P[ipart].Pos[2] - P[jpart].Pos[2];

                if ( Problem.Periodic ) {

                    while ( dx > boxhalf[0] ) { // closest image
                        dx -= boxsize[0];
                    }

                    while ( dx < -boxhalf[0] ) {
                        dx += boxsize[0];
                    }

                    while ( dy > boxhalf[1] ) {
                        dy -= boxsize[1];
                    }

                    while ( dy < -boxhalf[1] ) {
                        dy += boxsize[1];
                    }

                    while ( dz > boxhalf[2] ) {
                        dz -= boxsize[2];
                    }

                    while ( dz < -boxhalf[2] ) {
                        dz += boxsize[2];
                    }
                }

                float r2 = ( dx * dx + dy * dy + dz * dz );

                Assert ( r2 > 0,
                         "Found two particles %d & %d at the same location. "
                         "Consider increasing the space between your density field"
                         " and the box boundaries.", ipart, jpart );

                float h = 0.5 * ( hsml[ipart] + hsml[jpart] );

                if ( r2 > p2 ( h ) ) {
                    continue ;
                }

                float r = sqrt ( r2 );
                // * p3(h) is a for legacy reasons - at some point retune the code to work without it
                // norm_hsml also plays a minor role with that
#ifdef TWO_DIM
                double kernel_fac = p2 ( h );
#else
                double kernel_fac = p3 ( h );
#endif
                float wk = sph_kernel ( r, h ) * kernel_fac;

                // Also 1/3 for 2D since density contrast is not affected by dimensionality
                const double dens_contrast = pow ( SphP[ipart].Rho_Model / rho_mean, 1 / 3 );

                delta[0][ipart] += step / dens_contrast * hsml[ipart] * wk * dx / r;
                delta[1][ipart] += step / dens_contrast * hsml[ipart] * wk * dy / r;
#ifndef TWO_DIM
                delta[2][ipart] += step / dens_contrast * hsml[ipart] * wk * dz / r;
#endif
            }

        }

        int cnt = 0, cnt1 = 0, cnt2 = 0, cnt3 = 0;

        #pragma omp parallel for shared(delta,P) reduction(+:cnt,cnt1,cnt2,cnt3)
        for ( int ipart = 0; ipart < nPart; ipart++ ) { // move particles

            const float d = sqrt ( p2 ( delta[0][ipart] ) + p2 ( delta[1][ipart] ) + p2 ( delta[2][ipart] ) );

            const float h = SphP[ipart].Hsml;

#ifdef TWO_DIM
            float d_mps = pow ( pi * p2 ( h ) / DESNNGB, 1.0 / 3.0 );
#else
            float d_mps = pow ( fourpithird * p3 ( h ) / DESNNGB, 1.0 / 3.0 );
#endif

            if ( d > 1 * d_mps ) {
                ++cnt;
            }
            if ( d > 0.1 * d_mps ) {
                ++cnt1;
            }
            if ( d > 0.01 * d_mps ) {
                ++cnt2;
            }
            if ( d > 0.001 * d_mps ) {
                ++cnt3;
            }

            P[ipart].Pos[0] += delta[0][ipart]; // push !
            P[ipart].Pos[1] += delta[1][ipart];
            P[ipart].Pos[2] += delta[2][ipart];

            while ( P[ipart].Pos[0] < 0 ) { // keep it in the box
                P[ipart].Pos[0] += boxsize[0];
            }

            while ( P[ipart].Pos[0] > boxsize[0] ) {
                P[ipart].Pos[0] -= boxsize[0];
            }

            while ( P[ipart].Pos[1] < 0 ) {
                P[ipart].Pos[1] += boxsize[1];
            }

            while ( P[ipart].Pos[1] > boxsize[1] ) {
                P[ipart].Pos[1] -= boxsize[1];
            }

            while ( P[ipart].Pos[2] < 0 ) {
                P[ipart].Pos[2] += boxsize[2];
            }

            while ( P[ipart].Pos[2] > boxsize[2] ) {
                P[ipart].Pos[2] -= boxsize[2];
            }
        }

        double moveMps[4];
        moveMps[0] = cnt * 100. / Param.Npart;
        moveMps[1] = cnt * 100. / Param.Npart;
        moveMps[2] = cnt * 100. / Param.Npart;
        moveMps[3] = cnt * 100. / Param.Npart;

        printf ( "        Del %g%% > Dmps; %g%% > Dmps/10; %g%% > Dmps/100; %g%% > Dmps/1000\n",
                 moveMps[0], moveMps[1], moveMps[2], moveMps[3] );

#ifdef OUTPUT_DIAGNOSTICS
        struct Quadruplet errorQuad;
        errorQuad.min = errMin;
        errorQuad.max = errMax;
        errorQuad.mean = errMean;
        errorQuad.sigma = errSigma;

        const struct Quadruplet deltaQuad = calculateStatsOn ( delta, Param.Npart );

        writeIterationDiagnostics ( it, &errorQuad, errDiff, moveMps, &deltaQuad );
#endif

        if (   ( cnt * 100. / Param.Npart < Param.LimitMps[0] )
                || ( cnt1 * 100. / Param.Npart < Param.LimitMps[1] )
                || ( cnt2 * 100. / Param.Npart < Param.LimitMps[2] )
                || ( cnt3 * 100. / Param.Npart < Param.LimitMps[3] ) ) {
            break;
        }

        if ( cnt1 > last_cnt ) { // force convergence if distribution doesnt tighten
            step *= Param.StepReduction;
        }

        last_cnt = cnt1;
    }

    Free ( hsml );
    Free ( delta[0] );
    Free ( delta[1] );
    Free ( delta[2] );

    printf ( "done\n\n" );
    fflush ( stdout );

    return ;
}

void writeStepFile ( int it )
{

    char problem_name[CHARBUFSIZE] = "";
    char wvt_stepnumber[CHARBUFSIZE] = "";
    char wvt_stepname[CHARBUFSIZE] = "";
    sprintf ( problem_name, "%s", Problem.Name );

    strcpy ( wvt_stepname, problem_name );
    sprintf ( wvt_stepnumber, "_%03d", it );
    strcat ( wvt_stepname, wvt_stepnumber );
    sprintf ( Problem.Name, "%s", wvt_stepname );
    printf ( "Writing file %s\n", Problem.Name );
    Write_output ( 0 ); // not verbose
    sprintf ( Problem.Name, "%s", problem_name );
}
