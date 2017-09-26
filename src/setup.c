#include "globals.h"
#include "io.h"

/* Standard Gadget units */
const double ULength = 3.08568025e21; 	// kpc in cgs
const double UMass = 1.989e43;			// 10^10 Msol in cgs

float zeroFunction ( const int ipart );
void zeroFunctionVec ( const int ipart, float *out );

void setupProblemFromModel ( const int Flag, const int Subflag );

void mpartFromIntegral();

void Setup()
{
    size_t nBytes = Param.Npart * sizeof ( *P ); // allocate particles
    P = Malloc ( nBytes );
    memset ( P, 0, nBytes );

    nBytes = Param.Npart * sizeof ( *SphP );
    SphP = Malloc ( nBytes );
    memset ( SphP, 0, nBytes );

    setupProblemFromModel ( Param.Problem_Flag, Param.Problem_Subflag );
    mpartFromIntegral();

    printf ( "Problem %d.%d : %s \n"
             "   Npart: %d \n"
             "   Mpart: %g \n"
             "   Boxsize: %g x %g x %g \n"
             "   Periodic: %d \n\n",
             Param.Problem_Flag, Param.Problem_Subflag, Problem.Name,
             Param.Npart, Problem.Mpart, Problem.Boxsize[0],
             Problem.Boxsize[1], Problem.Boxsize[2], Problem.Periodic );
}

void setupProblemFromModel ( const int Flag, const int Subflag )
{
    Density_Func_Ptr = &zeroFunction;
    U_Func_Ptr = &zeroFunction;
    Velocity_Func_Ptr = &zeroFunctionVec;
    Magnetic_Field_Func_Ptr = &zeroFunctionVec;

    Problem.Mpart = 1; // required to renormalize later

    Problem.Periodic = true; // standard settings
    Problem.Rho_Max = 1.0;
    Problem.Boxsize[0] = Problem.Boxsize[1] = Problem.Boxsize[2] = 1;

    switch ( Flag ) {

    case 0:

        switch ( Subflag ) {

        case 0:
            setup_Constant_Density();
            break;

        case 1:
            setup_TopHat_Density();
            break;

        case 2:
            setup_Sawtooth_Density();
            break;

        case 3:
            setup_SineWave_Density();
            break;

        default:
            Assert ( false, "Effect %d.%d not implemented", Flag, Subflag );
            break;
        }

        break;

    case 1:

        Problem.Periodic = false;

        switch ( Subflag ) {

        case 0:
            setup_Gradient_Density();
            break;

        default:
            Assert ( false, "Effect %d.%d not implemented", Flag, Subflag );
            break;
        }

        break;

    case 2:

        switch ( Subflag ) {

        case 0:
            setup_Magneticum_Density();
            break;

        case 1:
            setup_Png_Density();
            break;

        case 2:
            setup_Grid_Densiy();
            break;

        case 3:
            Assert ( false, "Effect %d.%d not implemented", Flag, Subflag );
            break;

        default:
            Assert ( false, "Effect %d.%d not implemented", Flag, Subflag );
            break;
        }

        break;

    case 3:
        setup_Double_Shock ( Subflag );
        break;

    case 4:
        setup_User_Density ( Subflag );
        break;

    default:
        Assert ( false, "Effect %d.%d not implemented", Flag, Subflag );
        break;
    }

    Assert ( Problem.Boxsize[0] >= Problem.Boxsize[1] &&
             Problem.Boxsize[0] >= Problem.Boxsize[2],
             "Boxsize[0] has to be largest for ngb finding to work." );

    return ;
}

float zeroFunction ( const int ipart )
{
    return 0;
}

void zeroFunctionVec ( const int ipart, float *out )
{
    out[0] = out[1] = out[2] = 0;

    return;
}

void mpartFromIntegral()
{
    const int N = 1ULL << 9;
    const double dx = Problem.Boxsize[0] / N;
    const double dy = Problem.Boxsize[1] / N;
#ifndef TWO_DIM
    const double dz = Problem.Boxsize[2] / N;
#endif

    double tot_mass = 0.0;

    for ( int i = 0; i < N; i++ ) {
        P[0].Pos[0] = ( i + 0.5 ) * dx;

        for ( int j = 0; j < N; j++ ) {
            P[0].Pos[1] = ( j + 0.5 ) * dy;

#ifdef TWO_DIM
            P[0].Pos[2] = 0.0;
            tot_mass += Density_Func_Ptr ( 0 ) * dx * dy;
#else
            for ( int k = 0; k < N; k++ ) {
                P[0].Pos[2] = ( k + 0.5 ) * dz;
                tot_mass += Density_Func_Ptr ( 0 ) * dx * dy * dz;
            }
#endif
        }
    }

    Problem.Mpart = tot_mass / Param.Npart;

    Assert ( Problem.Mpart > 0,
             "Particle mass has to be finite, have %g", Problem.Mpart );

    return ;
}
