float betamodel ( const float r );
float GalaxyCluster_Density ( const int ipart );
float GalaxyCluster_Mass ( const int ipart );
void Setup_Gas_Mass_Profile();
float Gas_Mass_Profile ( const float r_in );
float Inverted_Gas_Mass_Profile ( float M );

void GalaxyCluster_Velocity ( const int ipart, float out[3] );
float GalaxyCluster_U ( const int ipart );

extern struct HaloProperties {
    float Rho0;
    float Beta;
    float Rcore;
    float R_Sample;
} Halo;
