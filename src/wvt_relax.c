#include "globals.h"
#include "tree.h"

#define WVTNNGB DESNNGB // 145 for WC2 that equals WC6

#define TREEBUILDFREQUENCY 1

#define MAXITER 64
#define MINITER 8
#define ERRDIFF_LIMIT 0.005
#define ERRMEAN_LIMIT 0.05
#define ERRMAX_LIMIT 0.01

#define STEP_DIVISOR 5.0
#define STEP_DECLINE 1.0
#define GLOBAL_STEP_FAC_DECLINE 0.9

int Find_ngb_simple(const int ipart,  const float hsml, int *ngblist);
int ngblist[NGBMAX] = { 0 }, Ngbcnt ;

static inline float sph_kernel_M4(const float r, const float h);
static inline double sph_kernel_WC2(const float r, const float h);
static inline double sph_kernel_WC6(const float r, const float h);
static inline float gravity_kernel(const float r, const float h);

void writeStepFile(int it);

static int global_it = -1;
static double global_step_fac = 1.0;

/* Settle SPH particle with weighted Voronoi tesselations (Diehl+ 2012).
 * Here hsml is not the SPH smoothing length, but is related to a local
 * metric defined ultimately by the density model.
 * Relaxation is done in units of the boxsize, hence the box volume is 1
 * Return code true means that rerun could be usefull */

bool Regularise_sph_particles()
{
    bool returnCode = false;

    const int nPart = Param.Npart;

    const double boxsize[3] = { Problem.Boxsize[0], Problem.Boxsize[1],
								Problem.Boxsize[2]};
    const double boxhalf[3] = { boxsize[0]/2, boxsize[1]/2, boxsize[2]/2, };
	const double boxinv[3] = { 1/boxsize[0], 1/boxsize[1], 1/boxsize[2] };

    printf("Starting iterative SPH regularisation \n"
            "   max %d iterations, tree update every %d iterations\n"
            "   stop at errdiff < %g%%, errmean < %g%%, errmax < %g%%   \n\n",
            MAXITER, TREEBUILDFREQUENCY, ERRDIFF_LIMIT*100, ERRMEAN_LIMIT*100,
			ERRMAX_LIMIT*100); fflush(stdout);

    float *hsml = NULL;
    size_t nBytes = nPart * sizeof(*hsml);
    hsml = Malloc(nBytes);

    float *delta[3] = { NULL };
    nBytes = nPart * sizeof(**delta);
    delta[0] = Malloc(nBytes);
    delta[1] = Malloc(nBytes);
    delta[2] = Malloc(nBytes);

    int it = global_it;

	double npart_1D = pow(nPart, 1.0/3.0);

    double step[3] = {  Problem.Boxsize[0] / npart_1D / STEP_DIVISOR * global_step_fac,
                        Problem.Boxsize[1] / npart_1D / STEP_DIVISOR * global_step_fac,
                        Problem.Boxsize[2] / npart_1D / STEP_DIVISOR * global_step_fac } ;
#ifdef SPH_CUBIC_SPLINE
	step[0] *= 6;
	step[1] *= 6;
	step[2] *= 6;
#endif // SPH_CUBIC_SPLINE

    double errLast = DBL_MAX, errLastTree = DBL_MAX;
    double errDiff = DBL_MAX, errDiffLast = DBL_MAX;

    const double volume = Problem.Boxsize[0] * Problem.Boxsize[1] * Problem.Boxsize[2];
    const double rho_mean = nPart * Problem.Mpart / volume;

    for (;;) {

        if ((it++ % TREEBUILDFREQUENCY) == 0)
            Find_sph_quantities();

#ifdef SAVE_WVT_STEPS
        writeStepFile(it);
#endif

        if (it-global_it+1 > MAXITER) {

            printf("Reached max iterations - ");
            returnCode = true;
            break;
        }

        int nIn = 0;
        double  errMax = 0, errMean = 0;

		#pragma omp parallel for reduction(+:errMean,nIn) reduction(max:errMax)
        for (int  ipart = 0; ipart < nPart; ipart++) { // get error

            float rho = (*Density_Func_Ptr) (ipart);

            float err = fabs(SphP[ipart].Rho-rho) / rho;

            errMax = fmax(err, errMax);

            errMean += err;

            nIn++;
        }

        errMean /= nIn;

        errDiff = (errLast - errMean) / errMean;

        printf("   #%02d: Err max=%3g mean=%03g diff=%03g"
                " step=%g %g %g \n", it, errMax, errMean,errDiff, 
				step[0], step[1], step[2]);

        if (fabs(errDiff) < ERRDIFF_LIMIT && it-global_it+1 > MINITER) {

            printf("Mean error hardly improves anymore - ");
            returnCode = true;
            break;
        }

        if ((fabs(errMean) < ERRMEAN_LIMIT &&
			fabs(errMax) < ERRMAX_LIMIT) && it-global_it+1 > MINITER) {

            printf("Achieved desired error criterion - ");
            returnCode = false;
            break;
        }

        if ((errDiff < 0) && (errDiffLast < 0) && (it-global_it+1 > MINITER)) { //stop if worse

            printf("Convergence flipped - ");
            returnCode = true;//false;
            break;
        }

        if ((fabs(errDiff) < 0.01)) { // force convergence
            step[0] *= STEP_DECLINE;
            step[1] *= STEP_DECLINE;
            step[2] *= STEP_DECLINE;
		}

        errLast = errMean;
        errDiffLast = errDiff;

        double vSphSum = 0; // total volume defined by hsml

		#pragma omp parallel for shared(hsml) reduction(+:vSphSum)
        for (int ipart = 0; ipart < nPart; ipart++) { // find hsml

            float rho = (*Density_Func_Ptr) (ipart);

            SphP[ipart].Rho_Model= rho;

            hsml[ipart] = pow(WVTNNGB*Problem.Mpart/rho/fourpithird, 1./3.);

            vSphSum += p3(hsml[ipart]);
        }

        float norm_hsml = pow(WVTNNGB/vSphSum/fourpithird , 1.0/3.0);

		#pragma omp parallel for
        for (int ipart = 0; ipart < nPart; ipart++)
            hsml[ipart] *= norm_hsml;

		#pragma omp parallel for shared(delta, hsml, P) \
        	schedule(dynamic, nPart/Omp.NThreads/256)
        for (int ipart = 0; ipart < nPart; ipart++) {

            delta[0][ipart] = delta[1][ipart] = delta[2][ipart] = 0;

            int ngblist[NGBMAX] = { 0 };
            int ngbcnt = Find_ngb_tree(ipart, hsml[ipart], ngblist);
            //int ngbcnt = Find_ngb_simple(ipart, hsml[ipart], ngblist);

            const float rho = (*Density_Func_Ptr) (ipart);
            const float mean_dist = pow(Problem.Mpart / rho / DESNNGB, 1.0/3.0);
            const float d_max = 1.0 * mean_dist;
            for (int i = 0; i < ngbcnt; i++) { // neighbour loop

                int jpart = ngblist[i];

                if (ipart == jpart)
                    continue;

                float dx = P[ipart].Pos[0] - P[jpart].Pos[0];
                float dy = P[ipart].Pos[1] - P[jpart].Pos[1];
                float dz = P[ipart].Pos[2] - P[jpart].Pos[2];

                if (Problem.Periodic) {

                    dx = dx > boxhalf[0] ? dx-boxsize[0] : dx; // closest image
                    dy = dy > boxhalf[1] ? dy-boxsize[1] : dy;
                    dz = dz > boxhalf[2] ? dz-boxsize[2] : dz;

                    dx = dx < -boxhalf[0] ? dx+boxsize[0] : dx;
                    dy = dy < -boxhalf[1] ? dy+boxsize[1] : dy;
                    dz = dz < -boxhalf[2] ? dz+boxsize[2] : dz;
                }

                float r2 = (dx*dx + dy*dy + dz*dz);

				Assert(r2 > 0, 
					"Found two particles %d & %d at the same location. "
					"Consider increasing the space between your density field"
					" and the box boundaries.", ipart, jpart);

                float h = 0.5 * (hsml[ipart] + hsml[jpart]);

                if (r2 > p2(h))
                    continue ;

                float r = sqrt(r2);

#ifdef SPH_CUBIC_SPLINE
                float wk = sph_kernel_WC2(r, h);
#else
                float wk = sph_kernel_WC6(r, h);
#endif

                const double dens_contrast = pow(SphP[ipart].Rho_Model/rho_mean, 1/3);

                delta[0][ipart] += step[0] / dens_contrast * hsml[ipart] * wk * dx / r;
                delta[1][ipart] += step[1] / dens_contrast * hsml[ipart] * wk * dy / r;
                delta[2][ipart] += step[2] / dens_contrast * hsml[ipart] * wk * dz / r;
            }
            delta[0][ipart] = fmin(delta[0][ipart], d_max);
            delta[1][ipart] = fmin(delta[1][ipart], d_max);
            delta[2][ipart] = fmin(delta[2][ipart], d_max);
        }

        int cnt = 0, cnt1 = 0, cnt2 = 0;

		#pragma omp parallel for shared(delta,P) \
        	reduction(+:cnt) reduction(+:cnt1) reduction(+:cnt2)
        for (int ipart = 0; ipart < nPart; ipart++) { // move particles

            float rho = (*Density_Func_Ptr) (ipart);

            float d = sqrt(p2(delta[0][ipart])
                    + p2( delta[1][ipart]) + p2( delta[2][ipart]));

            float d_mean = pow(Problem.Mpart / rho / DESNNGB, 1.0/3.0);

            if (d > 1 * d_mean)
                cnt++;
            if (d > 0.1 * d_mean)
                cnt1++;
            if (d > 0.01 * d_mean)
                cnt2++;

            P[ipart].Pos[0] += delta[0][ipart]; // push !
            P[ipart].Pos[1] += delta[1][ipart];
            P[ipart].Pos[2] += delta[2][ipart];

            while (P[ipart].Pos[0] < 0) // keep it in the box
                P[ipart].Pos[0] += boxsize[0];

            while (P[ipart].Pos[0] > boxsize[0])
                P[ipart].Pos[0] -= boxsize[0];

            while (P[ipart].Pos[1] < 0)
                P[ipart].Pos[1] += boxsize[1];

            while (P[ipart].Pos[1] > boxsize[1])
                P[ipart].Pos[1] -= boxsize[1];

            while (P[ipart].Pos[2] < 0)
                P[ipart].Pos[2] += boxsize[2];

            while (P[ipart].Pos[2] > boxsize[2])
                P[ipart].Pos[2] -= boxsize[2];
        }

        printf("        Del %g%% > dp; %g%% > dp/10; %g%% > dp/100\n", 
				cnt*100./Param.Npart, cnt1*100./Param.Npart, cnt2*100./Param.Npart);
    }

    global_it = it;
    global_step_fac *= GLOBAL_STEP_FAC_DECLINE;

    Free(hsml); Free(delta[0]); Free(delta[1]); Free(delta[2]);

    printf("done\n\n"); fflush(stdout);

    return returnCode;
}

void writeStepFile(int it) {

    char problem_name[CHARBUFSIZE] = "";
    char wvt_stepnumber[CHARBUFSIZE] = "";
    char wvt_stepname[CHARBUFSIZE] = "";
    sprintf(problem_name, Problem.Name);

    strcpy(wvt_stepname, problem_name);
    sprintf(wvt_stepnumber, "_%03d", it);
    strcat(wvt_stepname, wvt_stepnumber);
    sprintf(Problem.Name, wvt_stepname);
    printf("Writing file %s\n", Problem.Name);
    Write_output(0);  // not verbose
    sprintf(Problem.Name, problem_name);
}

static inline double sph_kernel_WC2(const float r, const float h)
{
    const float u= r/h;
    const float t = 1-u;

    return 21/(2*pi)*t*t*t*t*(1+4*u);
}

static inline float gravity_kernel(const float r, const float h)
{
    const float epsilon = 0.1;
    const float offset = h / (h + epsilon);
    const float val = h / (r + epsilon) - offset;

    return val * val;
}

static inline double sph_kernel_WC6(const float r, const float h)
{
    const double u = r/h;
    const double t = 1-u;

    return 1365.0/(64*pi) *t*t*t*t*t*t*t*t*(1+8*u + 25*u*u + 32*u*u*u);
}

static inline float sph_kernel_M4(const float r, const float h) // cubic spline
{
    double wk = 0;
    double u = r/h;

    if(u < 0.5)
        wk = (2.546479089470 + 15.278874536822 * (u - 1) * u * u);
    else
        wk = 5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u);

    return wk/p3(h);
}

