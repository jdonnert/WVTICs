#include "../globals.h"
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_vector.h>


void setup_Sedov_Blast()
{
    Problem.Boxsize[0] = 3;
    Problem.Boxsize[1] = 3;
    Problem.Boxsize[2] = 3;

    sprintf ( Problem.Name, "IC_SedovBlast" );

    const double rho = 1.24E7;

    Problem.Rho_Max = rho;

    Density_Func_Ptr = &Sedov_Blast_Density;
    U_Func_Ptr = &Sedov_Blast_U;
    PostProcessing_Func_Ptr = &Sedov_Blast_PostProcessing;
}

float Sedov_Blast_Density ( const int ipart , const double bias )
{
    return 1.24E7;
}

// This function calculates the distance from each particle to the center of the Box. It returns the distance to the NNpart times furthest neighbour of zero.

float Sedov_Blast_abs ()
{

    int NNpart = 32;

    gsl_vector *abs_of_zero = gsl_vector_alloc ( Param.Npart );

    for ( int i = 0; i < Param.Npart; i++ ) {

        double r = sqrt ( ( P[i].Pos[0] - 0.5 * Problem.Boxsize[0] ) * ( P[i].Pos[0] - 0.5 * Problem.Boxsize[0] ) + ( P[i].Pos[1] * 0.5 * Problem.Boxsize[1] ) * ( P[i].Pos[1] - 0.5 * Problem.Boxsize[1] )
                          + ( P[i].Pos[2] - 0.5 * Problem.Boxsize[2] ) * ( P[i].Pos[2] - 0.5 * Problem.Boxsize[2] ) );
        gsl_vector_set ( abs_of_zero, i, r );
    }

    gsl_sort_vector ( abs_of_zero );

    double dist = gsl_vector_get ( abs_of_zero, NNpart - 1 );

    gsl_vector_free ( abs_of_zero );

    return dist;

}

#ifdef KINETIC_SEDOV
float Sedov_Blast_kinetic ()
{

    int NNpart = 32;
    double Radius = partpos[3][32];

    float maxDistance = Sedov_Blast_abs ();

    if ( Radius <= maxDistance ) {

        for ( int i =; i < 3; i++ ) {

            for ( int j = 0; j < NNpart, j++ ) {

                double pos = P[i].Pos[j];
                partpos[i][j] = pos;

            }

        }

    }

}

#endif // KINETIC_SEDOV


float Sedov_Blast_U ( const int ipart)
{
    return 0.0;
}


//! @todo improvement: use gsl_sort_vector_index
void Sedov_Blast_PostProcessing ()
{
    /*const double u_sn = 4.18971E5;
    int sn_count = pow ( Param.Npart / 3200., 3.0 );
    sn_count = min ( sn_count, 1 );
    //! @todo assign u_sn to sn_count innermost particles*/

    float maxDistance = Sedov_Blast_abs ();

    for ( int ipart = 0; ipart < Param.Npart; ++ipart ) {


        const double x = P[ipart].Pos[0] - Problem.Boxsize[0] * 0.5;
        const double y = P[ipart].Pos[1] - Problem.Boxsize[1] * 0.5;
        const double z = P[ipart].Pos[2] - Problem.Boxsize[2] * 0.5;
        double Radius = sqrt ( x * x + y * y + z * z );
        // int NNpart = 32;
        const double  Supernova_Energy = 0.00502765; // This is the value of one supernova energy per unit mass in Gadget units 1e51 erg
        // const double hsml = SphP[i].Hsml;

        if ( Radius <=  maxDistance ) {
            SphP[ipart].U = Supernova_Energy;
        }
    }

}
