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

float GalaxyCluster_Density ( const int ipart , const double bias )
{
    const float x = P[ipart].Pos[0] - 0.5 * Problem.Boxsize[0];
    const float y = P[ipart].Pos[1] - 0.5 * Problem.Boxsize[1];;
    const float z = P[ipart].Pos[2] - 0.5 * Problem.Boxsize[2];;

    const float r = sqrt ( x * x + y * y + z * z );

    return betamodel ( r );
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
