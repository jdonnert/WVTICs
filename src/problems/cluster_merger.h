float betamodel_merger( float r, int cluster );
float ClusterMerger_Density ( const int ipart , const double bias );

void ClusterMerger_Velocity ( const int ipart, float out[3] );
float ClusterMerger_U ( const int ipart );

extern struct HaloProperties {
    float Pos[3];
    float Rho0;
    float Beta;
    float Rcore;
} Halo;

typedef struct {
    int Cluster;
} ParentCluster;
