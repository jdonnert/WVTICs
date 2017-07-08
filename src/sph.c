#include "globals.h"
#include "sph.h"
#include "tree.h"
#include "kernel.h"

#ifdef TWO_DIM
#define HSML_FACTOR 1.15
#else
#define HSML_FACTOR 1.24
#endif

extern void Find_sph_quantities()
{
    Sort_Particles_By_Peano_Key();

    Build_Tree();

    #pragma omp parallel for shared(SphP, P) \
    schedule(dynamic, Param.Npart/Omp.NThreads/64)
    for ( int ipart = 0; ipart < Param.Npart; ipart++ ) {

        float hsml = SphP[ipart].Hsml;

        if ( hsml == 0 ) {
            hsml = 2 * Guess_hsml ( ipart, DESNNGB );    // always too large
        }

        Assert ( isfinite ( hsml ), "hsml not finite ipart=%d parent=%d \n",
                 ipart, P[ipart].Tree_Parent );

        float dRhodHsml = 0;
        float rho = 0;

        for ( ;; ) {

            int ngblist[NGBMAX] = { 0 };
            int ngbcnt = Find_ngb ( ipart, hsml, ngblist );

            if ( ngbcnt == NGBMAX ) { // prevent overflow of ngblist

                hsml /= HSML_FACTOR;

                continue;
            }

            if ( ngbcnt < DESNNGB ) {

                hsml *= HSML_FACTOR;

                continue;
            }

            bool part_done = Find_hsml ( ipart, ngblist, ngbcnt, &dRhodHsml,
                                         &hsml, &rho );

            if ( ngbcnt < DESNNGB && ( !part_done ) ) {
                hsml *= HSML_FACTOR;
            }

            if ( part_done ) {
                break;
            }
        }

        float varHsmlFac = 1.0 / ( 1 + hsml / ( 3 * rho ) * dRhodHsml );

        SphP[ipart].Hsml = hsml;
        SphP[ipart].Rho = rho;
        SphP[ipart].VarHsmlFac = varHsmlFac;
    }

    return;
}

/* solve SPH continuity eq via Newton-Raphson, bisection and tree search */

extern bool Find_hsml ( const int ipart, const int *ngblist, const int ngbcnt,
                        float *dRhodHsml_out, float *hsml_out, float *rho_out )
{
    const double boxsize[3] = { Problem.Boxsize[0], Problem.Boxsize[1],
                                Problem.Boxsize[2]
                              };
    const double boxhalf[3] = { boxsize[0] / 2, boxsize[1] / 2, boxsize[2] / 2, };

    double upper = *hsml_out * sqrt3;
    double lower = 0;

    double hsml = *hsml_out; //lower + 0.5*(upper-lower);
    double rho = 0, dRhodHsml = 0;

    int it = 0;

    bool part_done = 0;

    for ( ;; ) {

        const double pos_i[3] = { P[ipart].Pos[0], P[ipart].Pos[1],
                                  P[ipart].Pos[2]
                                };

        double wkNgb = 0; // kernel weight number of neighbours

        rho = dRhodHsml = 0;

        it++;

        for ( int i = 0; i < ngbcnt; i++ ) {

            int jpart = ngblist[i];

            double dx = pos_i[0] - P[jpart].Pos[0];
            double dy = pos_i[1] - P[jpart].Pos[1];
            double dz = pos_i[2] - P[jpart].Pos[2];

            if ( Problem.Periodic ) {

                while ( dx > boxhalf[0] ) { // find closest image
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

            double r2 = dx * dx + dy * dy + dz * dz;

            if ( r2 > p2 ( hsml ) ) {
                continue ;
            }

            double r = sqrt ( r2 );

            double wk = sph_kernel ( r, hsml );
            double dwk = sph_kernel_derivative ( r, hsml );

#ifdef TWO_DIM
            wkNgb += pi * wk * p2 ( hsml );
#else
            wkNgb += fourpithird * wk * p3 ( hsml );
#endif

            rho += Problem.Mpart * wk;

            dRhodHsml += -Problem.Mpart * ( 3 / hsml * wk + r / hsml * dwk );
        }

        if ( it > 128 ) { // not enough neighbours ? -> hard exit
            break;
        }

        double ngbDev = fabs ( wkNgb - DESNNGB );

        if ( ngbDev < NNGBDEV ) {

            part_done = true;

            break;
        }

        if ( fabs ( upper - lower ) < 1e-4 ) { // find more neighbours !

            hsml *= 1.26; // double volume

            break;
        }

        if ( ngbDev < 0.5 * DESNNGB ) { // Newton Raphson

            double omega =  ( 1 + dRhodHsml * hsml / ( 3 * rho ) );

            double fac = 1 - ( wkNgb - DESNNGB ) / ( 3 * wkNgb * omega );

            fac = fmin ( HSML_FACTOR, fac ); // handle overshoot
            fac = fmax ( 1 / HSML_FACTOR, fac );

            hsml *= fac;

        } else {  // bisection

            if ( wkNgb > DESNNGB ) {
                upper = hsml;
            }

            if ( wkNgb < DESNNGB ) {
                lower = hsml;
            }

            hsml = pow ( 0.5 * ( p3 ( lower ) + p3 ( upper ) ), 1.0 / 3.0 );
        }
    } // for(;;)

    *hsml_out = ( float ) hsml;
    *rho_out = ( float ) rho;

    if ( part_done ) {
        *dRhodHsml_out = ( float ) dRhodHsml;

        double bias_corr = bias_correction ( hsml );
        *rho_out += bias_corr;
    }

    return part_done;
}
