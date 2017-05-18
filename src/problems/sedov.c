#include "../globals.h"

void setup_Sedov_Blast()
{
    Problem.Boxsize[0] = 6;
    Problem.Boxsize[1] = 6;
    Problem.Boxsize[2] = 6;

    sprintf ( Problem.Name, "IC_SedovBlast" );

    const double rho = 1.24E7;

    Problem.Rho_Max = rho;
    Problem.Mpart = rho * ( Problem.Boxsize[0] * Problem.Boxsize[1] * Problem.Boxsize[2] ) / Param.Npart;

    Density_Func_Ptr = &Sod_Shock_Density;
    U_Func_Ptr = &Sod_Shock_U;
    PostProcessing_Func_Ptr = &Sedov_Blast_PostProcessing;
}

float Sedov_Blast_Density ( const int ipart )
{
    return 1.24E7;
}

float Sedov_Blast_U ( const int ipart )
{
    return 0.0;
}

void Sedov_Blast_PostProcessing ()
{
    const double u_sn = 4.18971E5;
    int sn_count = pow ( Param.Npart / 3200., 3.0 );
    sn_count = min ( sn_count, 1 );
    //! @todo assign u_sn to sn_count innermost particles
}
