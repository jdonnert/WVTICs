#ifndef GLOBALS_H
#define GLOBALS_H

#define PROG_NAME "WVT Box"
#define VERSION "1.0"

/* C std lib */
#include <stdlib.h>             // system
#include <stdio.h>
#include <stdint.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <complex.h>

/* GNU Scientifc Library */
#include <gsl/gsl_math.h>
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_rng.h>

#include <omp.h>

/* Global Prototypes */
#include "macro.h"
#include "peano.h"
#include "proto.h"

/* Code parameters */

#define CHARBUFSIZE 512        // For any char buffer !
#define MAXTAGS 300            // In parameter file

#ifdef SPH_CUBIC_SPLINE

#define DESNNGB 50         // SPH kernel weighted number of neighbours
#define NNGBDEV 0.05       // error tolerance in SPH kernel weighted neighb.
#define NGBMAX (DESNNGB*16)  // size of neighbour list

#else

#define DESNNGB 295        // SPH kernel weighted number of neighbours
#define NNGBDEV 0.05       // error tolerance in SPH kernel weighted neighb.
#define NGBMAX (DESNNGB*8)  // size of neighbour list

#endif // SPH_CUBIC_SPLINE

/* mathematical constants */
#define pi             M_PI
#define sqrt2        M_SQRT2
#define sqrt3       1.73205080756887719
#define fourpithird 4.18879032135009765

extern struct OpenMP_infos{
    int NThreads;          // Number of openMP threads
    int ThreadID;          // Thread ID of this thread
    unsigned short Seed[3]; // random number seed: erand48(Omp.Seed)
} Omp;
#pragma omp threadprivate(Omp)

extern struct Parameters {
    int Npart;
	int Problem_Flag;
	int Problem_Subflag;
} Param;

extern struct ProblemParameters {
    char Name[CHARBUFSIZE];
	double Mpart;
    double Boxsize[3]; // [0] is ALWAYS the largest one !
	double Rho_Max;
	bool Periodic;
} Problem;

extern struct ParticleData{
    float Pos[3];
    float Vel[3];
    int32_t ID;
    int Type;
    peanoKey Key;
    int Tree_Parent;
} *P;

extern struct GasParticleData {
    float U;
    float Rho;
    float Hsml;
    float VarHsmlFac;
    float ID;
    float Rho_Model;
    float Bfld[3];
} *SphP;

extern float (*Density_Func_Ptr) (const int);
extern float (*U_Func_Ptr) (const int);
extern void (*Velocity_Func_Ptr) (const int, float *);
extern void (*Magnetic_Field_Func_Ptr) (const int, float *);

double G; // gravitational constant in code units

#endif
