#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include "../globals.h"
#include "galaxy_cluster.h"


#define NTABLE 256
#define GSL_SPLINE gsl_interp_cspline

#define Zero 1
#define Infinity 1e25

struct HaloProperties Halo = { 0 };

void setup_GalaxyCluster()
{
    Problem.Boxsize[0] = Problem.Boxsize[1] = Problem.Boxsize[2] = 1000;  // kpc
    Halo.Rho0 = 1e-26;  // g/cm^3
    Halo.Beta = 2.0 / 3.0;
    Halo.Rcore = 20;    // kpc
    Halo.R_Sample = Problem.Boxsize[0];

    sprintf ( Problem.Name, "IC_GalaxyCluster" );

    Problem.Rho_Max = Halo.Rho0;

    Setup_Gas_Mass_Profile();

    Density_Func_Ptr = &GalaxyCluster_Density;
    U_Func_Ptr = &GalaxyCluster_U;
    Velocity_Func_Ptr = &GalaxyCluster_Velocity;
}


float betamodel ( const float r )
{
    float rho0 = Halo.Rho0;
    float beta = Halo.Beta;
    float rc = Halo.Rcore;

    return rho0 * pow ( 1 + p2 ( r / rc ), -1.5 * beta );
}

float GalaxyCluster_Density ( const int ipart )
{
    const float x = P[ipart].Pos[0] - 0.5 * Problem.Boxsize[0];
    const float y = P[ipart].Pos[1] - 0.5 * Problem.Boxsize[1];;
    const float z = P[ipart].Pos[2] - 0.5 * Problem.Boxsize[2];;

    const float r = sqrt ( x * x + y * y + z * z );

    return betamodel ( r );
}


/* Cumulative Mass profile from the gas density profile by numerical
 * integration and cubic spline interpolation. This has to be called once
 * before the mass profile of a halo is used, to set the spline variables. */

static gsl_spline *M_Spline = NULL;
static gsl_interp_accel *M_Acc = NULL;
#pragma omp threadprivate(M_Spline, M_Acc)

static gsl_spline *Minv_Spline = NULL;
static gsl_interp_accel *Minv_Acc = NULL;
#pragma omp threadprivate(Minv_Spline, Minv_Acc)

float Gas_Mass_Profile ( const float r_in )
{
    double r = fmin ( r_in, Halo.R_Sample );

    return  gsl_spline_eval ( M_Spline, r, M_Acc );
}

float Inverted_Gas_Mass_Profile ( float M )
{
    return gsl_spline_eval ( Minv_Spline, M, Minv_Acc );
}

double m_integrant ( double r, void *param )
{
    return 4 * pi * r * r * betamodel ( r );
}

void Setup_Gas_Mass_Profile()
{
    double m_table[NTABLE] = { 0 };
    double r_table[NTABLE] = { 0 };
    double dr_table[NTABLE] = { 0 };

    double rmin = Zero;

    double rmax = Halo.R_Sample * 1.2; // include R_Sample

    double log_dr = ( log10 ( rmax / rmin ) ) / ( NTABLE - 1 );

    gsl_function gsl_F = { 0 };

    gsl_integration_workspace *gsl_workspace = NULL;
    gsl_workspace = gsl_integration_workspace_alloc ( NTABLE );

    for ( int i = 1; i < NTABLE; i++ ) {

        double error = 0;

        r_table[i] = rmin * pow ( 10, log_dr * i );

        gsl_F.function = &m_integrant;

        gsl_integration_qag ( &gsl_F, 0, r_table[i], 0, 1e-5, NTABLE,
                              GSL_INTEG_GAUSS61, gsl_workspace, &m_table[i], &error );

        if ( m_table[i] < m_table[i - 1] ) {
            m_table[i] = m_table[i - 1];    // integrator may fluctuate
        }

        //printf("%g \n", m_table[i]);
    }

    m_table[0] = 0;

    #pragma omp parallel
    {

        M_Acc  = gsl_interp_accel_alloc();

        M_Spline = gsl_spline_alloc ( GSL_SPLINE, NTABLE );
        gsl_spline_init ( M_Spline, r_table, m_table, NTABLE );

        Minv_Acc  = gsl_interp_accel_alloc();

        Minv_Spline = gsl_spline_alloc ( GSL_SPLINE, NTABLE );
        gsl_spline_init ( Minv_Spline, m_table, r_table, NTABLE );

    } // omp parallel

    return ;
}


float GalaxyCluster_Mass ( const int ipart )
{
    return 0.0;
}


void GalaxyCluster_Velocity ( const int ipart, float out[3] )
{
    out[0] = 0.0;
    out[1] = 0.0;
    out[2] = 0.0;
}

float GalaxyCluster_U ( const int ipart )
{
    return 0;
}
