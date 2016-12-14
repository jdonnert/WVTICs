#include "globals.h"
#include "tree.h"

#define WVTNNGB DESNNGB // 145 for WC2 that equals WC6

#define TREEBUILDFREQUENCY 1
#define NUMITER 64
#define ERRDIFF_LIMIT 0.005

int Find_ngb_simple(const int ipart,  const float hsml, int *ngblist);
int ngblist[NGBMAX] = { 0 }, Ngbcnt ;

static inline float sph_kernel_M4(const float r, const float h);
static inline double sph_kernel_WC2(const float r, const float h);
static inline double sph_kernel_WC6(const float r, const float h);
static inline float gravity_kernel(const float r, const float h);

/* Settle SPH particle with weighted Voronoi tesselations (Diehl+ 2012).
 * Here hsml is not the SPH smoothing length, but is related to a local
 * metric defined ultimately by the density model.
 * Relaxation is done in units of the boxsize, hence the box volume is 1 */

void Regularise_sph_particles()
{
    const int nPart = Param.Npart;

    const double boxsize[3] = { Problem.Boxsize[0], Problem.Boxsize[1],
								Problem.Boxsize[2]};
    const double boxhalf[3] = { boxsize[0]/2, boxsize[1]/2, boxsize[2]/2, };
	const double boxinv[3] = { 1/boxsize[0], 1/boxsize[1], 1/boxsize[2] };

#ifdef SAVE_WVT_STEPS
	char problem_name[CHARBUFSIZE] = "";
	char wvt_stepnumber[CHARBUFSIZE] = "";
	char wvt_stepname[CHARBUFSIZE] = "";
	sprintf(problem_name, Problem.Name);
#endif  // SAVE_WVT_STEPS

    printf("Starting iterative SPH regularisation \n"
            "   max %d iterations, tree update every %d iterations\n"
            "   stop at  errmax < %g%%   \n\n",
            NUMITER, TREEBUILDFREQUENCY, ERRDIFF_LIMIT*100); fflush(stdout);

    float *hsml = NULL;
    size_t nBytes = nPart * sizeof(*hsml);
    hsml = Malloc(nBytes);

    float *delta[3] = { NULL };
    nBytes = nPart * sizeof(**delta);
    delta[0] = Malloc(nBytes);
    delta[1] = Malloc(nBytes);
    delta[2] = Malloc(nBytes);

    int it = -1;

	double volume = Problem.Boxsize[0]*Problem.Boxsize[1]*Problem.Boxsize[2];
	double mean_part_sep = pow( volume / nPart, 1.0/3.0);
#ifdef SPH_CUBIC_SPLINE
    double step = 6 * mean_part_sep / 2;
#else
    double step = mean_part_sep / 2;
#endif // SPH_CUBIC_SPLINE

    double errLast = DBL_MAX, errLastTree = DBL_MAX;
    double errDiff = DBL_MAX, errDiffLast = DBL_MAX;

    for (;;) {

        if (it++ >= NUMITER) {
            printf("Reached max iterations - ");
            break;
        }

        if ((it % TREEBUILDFREQUENCY) == 0)
            Find_sph_quantities();

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
                " step=%g\n", it, errMax, errMean,errDiff, step);

#ifdef SAVE_WVT_STEPS
        strcpy(wvt_stepname, problem_name);
        sprintf(wvt_stepnumber, "_%03d", it);
        strcat(wvt_stepname, wvt_stepnumber);
        sprintf(Problem.Name, wvt_stepname);
        Write_output(0);  // not verbose
        sprintf(Problem.Name, problem_name);
#endif // SAVE_WVT_STEPS

        if (fabs(errDiff) < ERRDIFF_LIMIT && it > 32) { // at least iterate N times
            printf("Achieved desired error criterion - ");
            break;
        }

        if ((errDiff < 0) && (errDiffLast < 0) && (it > 10)) { // stop if worse
            printf("Convergence flipped - ");
            break;
        }

        if ((errDiff < 0.01) && (it > 1)) // force convergence
            step *= 0.8;

        errLast = errMean;
        errDiffLast = errDiff;

        double vSphSum = 0; // total volume defined by hsml

		#pragma omp parallel for shared(hsml) reduction(+:vSphSum)
        for (int ipart = 0; ipart < nPart; ipart++) { // find hsml

            float rho = (*Density_Func_Ptr) (ipart);

            SphP[ipart].Rho_Model= rho;

            hsml[ipart] = pow(WVTNNGB * Problem.Mpart /rho /fourpithird, 1./3.);

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

            int ngbcnt = Find_ngb_tree(ipart, hsml[ipart]*boxsize[0], ngblist);

            for (int i = 0; i < ngbcnt; i++) { // neighbour loop

                int jpart = ngblist[i];

                if (ipart == jpart)
                    continue;

                float dx = (P[ipart].Pos[0] - P[jpart].Pos[0]) * boxinv[0];
                float dy = (P[ipart].Pos[1] - P[jpart].Pos[1]) * boxinv[1];
                float dz = (P[ipart].Pos[2] - P[jpart].Pos[2]) * boxinv[2];

                dx = dx > 0.5 ? dx-1 : dx; // find closest image
                dy = dy > 0.5 ? dy-1 : dy;
                dz = dz > 0.5 ? dz-1 : dz;

                dx = dx < -0.5 ? dx+1 : dx;
                dy = dy < -0.5 ? dy+1 : dy;
                dz = dz < -0.5 ? dz+1 : dz;

                float r2 = (dx*dx + dy*dy + dz*dz);

                float h = 0.5 * (hsml[ipart] + hsml[jpart]);

                if (r2 > p2(h))
                    continue ;

                float r = sqrt(r2);
                float wk = sph_kernel_WC6(r, h);

                delta[0][ipart] += step * hsml[ipart] * wk * dx/r;
                delta[1][ipart] += step * hsml[ipart] * wk * dy/r;
                delta[2][ipart] += step * hsml[ipart] * wk * dz/r;
            }
        }

        int cnt = 0, cnt1 = 0, cnt2 = 0;

		#pragma omp parallel for shared(delta,P) \
        	reduction(+:cnt) reduction(+:cnt1) reduction(+:cnt2)
        for (int ipart = 0; ipart < nPart; ipart++) { // move particles

            float rho = (*Density_Func_Ptr) (ipart);

            float d = boxsize[0] * sqrt(p2(delta[0][ipart])
                    + p2( delta[1][ipart]) + p2( delta[2][ipart]));

            float meanPartSep = pow(Problem.Mpart / rho / DESNNGB, 1.0/3.0);

            if (d > 1 * meanPartSep)
                cnt++;
            if (d > 0.1 * meanPartSep)
                cnt1++;
            if (d > 0.01 * meanPartSep)
                cnt2++;

            P[ipart].Pos[0] += (float) (delta[0][ipart] * boxsize[0]); // push !
            P[ipart].Pos[1] += (float) (delta[1][ipart] * boxsize[1]);
            P[ipart].Pos[2] += (float) (delta[2][ipart] * boxsize[2]);

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
    }

    Free(hsml); Free(delta[0]); Free(delta[1]); Free(delta[2]);

    printf("done\n\n"); fflush(stdout);

    return ;
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

int Find_ngb_simple(const int ipart,  const float hsml, int *ngblist)
{
    const double boxsize[3] = { Problem.Boxsize[0], Problem.Boxsize[1],
								Problem.Boxsize[2]};
    const double boxhalf[3] = { boxsize[0]/2, boxsize[1]/2, boxsize[2]/2 };

    memset(ngblist, 0, NGBMAX*sizeof(*ngblist));

    int ngbcnt = 0;

    for (int jpart = 0; jpart < Param.Npart; jpart++) {

        float dx = (P[ipart].Pos[0] - P[jpart].Pos[0]);
        float dy = (P[ipart].Pos[1] - P[jpart].Pos[1]);
        float dz = (P[ipart].Pos[2] - P[jpart].Pos[2]);

        if (dx > boxhalf[0])    // find closest image
            dx -= boxsize[0];

        if (dy > boxhalf[1])
            dy -= boxsize[1];

        if (dz > boxhalf[2])
            dz -= boxsize[2];

        if (dx < -boxhalf[0])
            dx += boxsize[0];

        if (dy < -boxhalf[1])
            dy += boxsize[1];

        if (dz < -boxhalf[2])
            dz += boxsize[2];

        float r2 = dx*dx + dy*dy + dz*dz;

        if (r2 < hsml*hsml)
            ngblist[ngbcnt++] = jpart;

        if (ngbcnt == NGBMAX) {

			printf("WARNING, ngbcnt == %d, increase NGBMAX ! ",
					ngbcnt);

			break;
		}
    }

    return ngbcnt ;
}

