#include "globals.h"
#include "tree.h"
#include "kernel.h"

#define WVTNNGB DESNNGB

void writeStepFile ( int it );

/* Settle SPH particle with weighted Voronoi tesselations (Diehl+ 2012).
 * Here hsml is not the SPH smoothing length, but is related to a local
 * metric defined ultimately by the density model.
 * Relaxation is done in units of the boxsize, hence the box volume is 1
 * Return code true means that rerun could be usefull */

void Regularise_sph_particles()
{
    const int maxiter = 256;
    const double mps_frac = 5; 		// move this fraction of the mean particle sep
    const double step_red = 0.95; 	// force convergence at this rate
    const double bin_limits[3] = { -1, -1, 1};

    const int nPart = Param.Npart;

    const double boxsize[3] = { Problem.Boxsize[0], Problem.Boxsize[1],
                                Problem.Boxsize[2]
                              };
    const double boxhalf[3] = { boxsize[0] / 2, boxsize[1] / 2, boxsize[2] / 2, };

    const double median_boxsize = fmax ( boxsize[1], boxsize[2] ); // boxsize[0] is largest

    printf ( "Starting iterative SPH regularisation \n"
             "   Maxiter=%d, mps_frac=%g step_red=%g bin_limits=(%g,%g,%g)\n\n",
             maxiter, mps_frac, step_red, bin_limits[0], bin_limits[1], bin_limits[2] );
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

    double npart_1D = pow ( nPart, 1.0 / 3.0 );

    double step[3] = {  Problem.Boxsize[0] / npart_1D / mps_frac,
                        Problem.Boxsize[1] / npart_1D / mps_frac,
                        Problem.Boxsize[2] / npart_1D / mps_frac
                     } ;

    double errLast = DBL_MAX;
    double errDiff = DBL_MAX;

    const double volume = Problem.Boxsize[0] * Problem.Boxsize[1] * Problem.Boxsize[2];
    const double rho_mean = nPart * Problem.Mpart / volume;

    double last_cnt = DBL_MAX;

    for ( ;; ) {

        Find_sph_quantities();

        if ( it++ > maxiter ) {
            break;
        }

#ifdef SAVE_WVT_STEPS
        writeStepFile ( it );
#endif

        int nIn = 0;
        double  errMax = 0, errMean = 0;

        #pragma omp parallel for reduction(+:errMean,nIn) reduction(max:errMax)
        for ( int  ipart = 0; ipart < nPart; ipart++ ) { // get error

            float rho = ( *Density_Func_Ptr ) ( ipart );

            float err = fabs ( SphP[ipart].Rho - rho ) / rho;

            errMax = fmax ( err, errMax );

            errMean += err;

            nIn++;
        }

        errMean /= nIn;

        errDiff = ( errLast - errMean ) / errMean;

        printf ( "   #%02d: Err max=%3g mean=%03g diff=%03g"
                 " step=%g %g %g \n", it, errMax, errMean, errDiff,
                 step[0], step[1], step[2] );

        errLast = errMean;

        double vSphSum = 0; // total volume defined by hsml
        double max_hsml = 0;

        #pragma omp parallel for shared(hsml) reduction(+:vSphSum,max_hsml)
        for ( int ipart = 0; ipart < nPart; ipart++ ) { // find hsml

            float rho = ( *Density_Func_Ptr ) ( ipart );

            SphP[ipart].Rho_Model = rho;

            hsml[ipart] = pow ( WVTNNGB * Problem.Mpart / rho / fourpithird, 1. / 3. );

            vSphSum += p3 ( hsml[ipart] );

            max_hsml = max ( max_hsml, hsml[ipart] );
        }

        float norm_hsml = pow ( WVTNNGB / vSphSum / fourpithird , 1.0 / 3.0 ) * median_boxsize;

        #pragma omp parallel for
        for ( int ipart = 0; ipart < nPart; ipart++ ) {
            hsml[ipart] *= norm_hsml;
        }

        #pragma omp parallel for shared(delta, hsml, P) schedule(dynamic, nPart/Omp.NThreads/256)
        for ( int ipart = 0; ipart < nPart; ipart++ ) {

            delta[0][ipart] = delta[1][ipart] = delta[2][ipart] = 0;

            int ngblist[NGBMAX] = { 0 };
            int ngbcnt = Find_ngb_tree ( ipart, hsml[ipart], ngblist );
            //int ngbcnt = Find_ngb_simple(ipart, hsml[ipart], ngblist);

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
                float wk = sph_kernel ( r, h );

                const double dens_contrast = pow ( SphP[ipart].Rho_Model / rho_mean, 1 / 3 );

                delta[0][ipart] += step[0] / dens_contrast * hsml[ipart] / boxsize[0] * wk * dx / r ;
                delta[1][ipart] += step[1] / dens_contrast * hsml[ipart] / boxsize[1] * wk * dy / r ;
                delta[2][ipart] += step[2] / dens_contrast * hsml[ipart] / boxsize[2] * wk * dz / r ;
            }

        }

        int cnt = 0, cnt1 = 0, cnt2 = 0;

        #pragma omp parallel for shared(delta,P) reduction(+:cnt,cnt1,cnt2)
        for ( int ipart = 0; ipart < nPart; ipart++ ) { // move particles

            float rho = ( *Density_Func_Ptr ) ( ipart );

            float d = sqrt ( p2 ( delta[0][ipart] )
                             + p2 ( delta[1][ipart] ) + p2 ( delta[2][ipart] ) );

            float d_mps = pow ( Problem.Mpart / rho / DESNNGB, 1.0 / 3.0 );

            if ( d > 1 * d_mps ) {
                cnt++;
            }
            if ( d > 0.1 * d_mps ) {
                cnt1++;
            }
            if ( d > 0.01 * d_mps ) {
                cnt2++;
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

        printf ( "        Del %g%% > Dmps; %g%% > Dmps/10; %g%% > Dmps/100\n",
                 cnt * 100. / Param.Npart, cnt1 * 100. / Param.Npart, cnt2 * 100. / Param.Npart );

        if (   ( cnt * 100. / Param.Npart < bin_limits[0] )
                || ( cnt1 * 100. / Param.Npart < bin_limits[1] )
                || ( cnt2 * 100. / Param.Npart < bin_limits[2] ) ) {
            break;
        }

        if ( cnt1 > last_cnt ) { // force convergence if distribution doesnt tighten

            step[0] *= step_red;
            step[1] *= step_red;
            step[2] *= step_red;
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
    sprintf ( problem_name, Problem.Name );

    strcpy ( wvt_stepname, problem_name );
    sprintf ( wvt_stepnumber, "_%03d", it );
    strcat ( wvt_stepname, wvt_stepnumber );
    sprintf ( Problem.Name, wvt_stepname );
    printf ( "Writing file %s\n", Problem.Name );
    Write_output ( 0 ); // not verbose
    sprintf ( Problem.Name, problem_name );
}
