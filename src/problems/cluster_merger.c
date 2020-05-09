#include "../globals.h"
#include "cluster_merger.h"


struct HaloProperties Halo0 = { 0 };
struct HaloProperties Halo1 = { 0 };

static ParentCluster *Parent;

void setup_ClusterMerger()
{
    Problem.Boxsize[0] = 2000.0; // kpc
    Problem.Boxsize[1] = Problem.Boxsize[2] = 2000;  // kpc

    // Setup of cluster properties
    // -- Convention: Halo0 on the left!

    Halo0.Rho0   = 1e-26;  // g/cm^3
    Halo0.Beta   = 2.0 / 3.0;
    Halo0.Rcore  = 20;    // kpc
    Halo0.Pos[0] = 500.0;
    Halo0.Pos[1] = 1000.0;
    Halo0.Pos[2] = 1000.0;

    Halo1.Rho0  = 1e-26;  // g/cm^3
    Halo1.Beta  = 2.0 / 3.0;
    Halo1.Rcore = 15;    // kpc
    Halo1.Pos[0] = 1500.0;
    Halo1.Pos[1] = 1000.0;
    Halo1.Pos[2] = 1000.0;

    // allocate cluster parent array
    size_t nBytes = Param.Npart * sizeof ( ParentCluster );
    Parent = Malloc ( nBytes );
    memset ( Parent, 0, nBytes );

    sprintf ( Problem.Name, "IC_ClusterMerger" );

    // Cluster 0 needs to have the higher density!
    Problem.Rho_Max = Halo1.Rho0;

    Density_Func_Ptr = &ClusterMerger_Density;
    U_Func_Ptr = &ClusterMerger_U;
    Velocity_Func_Ptr = &ClusterMerger_Velocity;
}


float betamodel_merger ( float r, int cluster)
{
    float rho0, beta, rc;

    if ( cluster == 0)
    {
        rho0 = Halo0.Rho0;
        beta = Halo0.Beta;
        rc = Halo0.Rcore;
    }
    else
    {
        rho0 = Halo1.Rho0;
        beta = Halo1.Beta;
        rc = Halo1.Rcore;
    }


    return rho0 * pow ( 1 + p2 ( r / rc ), -1.5 * beta );
}

float ClusterMerger_Density ( const int ipart, const double bias )
{
    // Cluster 0
    float x = P[ipart].Pos[0] - Halo0.Pos[0];
    float y = P[ipart].Pos[1] - Halo0.Pos[1];
    float z = P[ipart].Pos[2] - Halo0.Pos[2];
    float r = sqrt ( x * x + y * y + z * z );
    float rho0 = betamodel_merger(r, 0);

    // Cluster 1
    x = P[ipart].Pos[0] - Halo1.Pos[0];
    y = P[ipart].Pos[1] - Halo1.Pos[1];
    z = P[ipart].Pos[2] - Halo1.Pos[2];
    r = sqrt ( x * x + y * y + z * z );
    float rho1 = betamodel_merger(r, 1);

    // define parent cluster via dominating density field
    if (rho0 > rho1)
        Parent[ipart].Cluster = 0;
    else
        Parent[ipart].Cluster = 1;

    return rho0 + rho1;
}

void ClusterMerger_Velocity ( const int ipart, float out[3] )
{
    // left cluster
    if ( !Parent[ipart].Cluster )
    {
        out[0] = 1000.0;
        out[1] = 0.0;
        out[2] = 0.0;
    }
    else // right cluster
    {
        out[0] = -2000.0;
        out[1] = 0.0;
        out[2] = 0.0;
    }

}

float ClusterMerger_U ( const int ipart )
{
    return 0;
}
