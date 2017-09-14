#include "globals.h"

/* Standard Gadget units */
const double ULength = 3.08568025e21; 	// kpc in cgs
const double UMass = 1.989e43;			// 10^10 Msol in cgs
const double UVel = 1e5;					// km/s in cgs

void zero_function_empty ();
float zero_function ( const int ipart );
void zero_function_vec ( const int ipart, float out[3] );

void setup_problem ( const int Flag, const int Subflag );
void mpart_from_integral();

void Setup()
{
    size_t nBytes = Param.Npart * sizeof ( *P ); // allocate particles
    P = Malloc ( nBytes );
    memset ( P, 0, nBytes );

    nBytes = Param.Npart * sizeof ( *SphP );
    SphP = Malloc ( nBytes );
    memset ( SphP, 0, nBytes );

    setup_problem ( Param.Problem_Flag, Param.Problem_Subflag );

    mpart_from_integral();

    printf ( "Problem %d.%d : %s \n"
             "   Npart: %d \n"
             "   Mpart: %g \n"
             "   Boxsize: %g x %g x %g \n"
             "   Periodic: %d \n\n",
             Param.Problem_Flag, Param.Problem_Subflag, Problem.Name,
             Param.Npart, Problem.Mpart, Problem.Boxsize[0],
             Problem.Boxsize[1], Problem.Boxsize[2], Problem.Periodic );
}

void setup_problem ( const int Flag, const int Subflag )
{
    Density_Func_Ptr = &zero_function;
    U_Func_Ptr = &zero_function;
    Velocity_Func_Ptr = &zero_function_vec;
    Magnetic_Field_Func_Ptr = &zero_function_vec;
    PostProcessing_Func_Ptr = &zero_function_empty;

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
        default:
            Assert ( false, "Effect %d.%d not implemented", Flag, Subflag );
            break;
        }
        break;

    case 3:
        setup_Double_Shock ( Subflag );
        break;

    case 4:
        switch ( Subflag ) {
        case 0:
            setup_Sod_Shock();
            break;
        case 1:
            setup_Sedov_Blast();
            break;
        case 2:
            setup_Kelvin_Helmholtz_Instability();
            break;
        case 4:
            setup_Blob();
            break;
        case 6:
            setup_Evrard_Collapse();
            break;
        case 7:
            setup_Zeldovich_Pancake();
            break;
        case 8:
            setup_Box();
            break;
        case 9:
            setup_Gresho_Vortex();
            break;
        case 10:
            setup_Exponential_Disk();
            break;
        case 11:
            setup_Boss();
            break;
        case 12:
            setup_GalaxyCluster();
            break;
        default:
            Assert ( false, "Effect %d.%d not implemented", Flag, Subflag );
            break;
        }
        break;

    case 5:
        switch ( Subflag ) {
        case 0:
            setup_Rotor();
            break;
        case 1:
            setup_Strong_Blast();
            break;
        case 2:
            setup_Orszag_Tang_Vortex();
            break;
        case 3:
            setup_Linear_Alfven_Wave();
            break;
        case 4:
            setup_Rayleigh_Taylor_Instability();
            break;
        case 5:
            set_Ryu_Jones_Parameters_1A();
            setup_Ryu_Jones_Shocktube();
            break;
        case 6:
            set_Ryu_Jones_Parameters_1B();
            setup_Ryu_Jones_Shocktube();
            break;
        case 7:
            set_Ryu_Jones_Parameters_2A();
            setup_Ryu_Jones_Shocktube();
            break;
        case 8:
            set_Ryu_Jones_Parameters_2B();
            setup_Ryu_Jones_Shocktube();
            break;
        case 9:
            set_Ryu_Jones_Parameters_3A();
            setup_Ryu_Jones_Shocktube();
            break;
        case 10:
            set_Ryu_Jones_Parameters_3B();
            setup_Ryu_Jones_Shocktube();
            break;
        case 11:
            set_Ryu_Jones_Parameters_4A();
            setup_Ryu_Jones_Shocktube();
            break;
        case 12:
            set_Ryu_Jones_Parameters_4B();
            setup_Ryu_Jones_Shocktube();
            break;
        case 13:
            set_Ryu_Jones_Parameters_4C();
            setup_Ryu_Jones_Shocktube();
            break;
        case 14:
            set_Ryu_Jones_Parameters_4D();
            setup_Ryu_Jones_Shocktube();
            break;
        case 15:
            set_Ryu_Jones_Parameters_5A();
            setup_Ryu_Jones_Shocktube();
            break;
        case 16:
            set_Ryu_Jones_Parameters_5B();
            setup_Ryu_Jones_Shocktube();
            break;
        default:
            Assert ( false, "MHD tests not implemented yet", Flag, Subflag );
            break;
        }
        break;

    case 6:
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

void zero_function_empty ()
{
    return;
}

float zero_function ( const int ipart )
{
    return 0;
}

void zero_function_vec ( const int ipart, float out[3] )
{
    out[0] = out[1] = out[2] = 0;

    return;
}

void mpart_from_integral()
{
    const int N = 1ULL << 9;
    const double dx = Problem.Boxsize[0] / N;
    const double dy = Problem.Boxsize[1] / N;
#ifndef TWO_DIM
    const double dz = Problem.Boxsize[2] / N;
#endif
    printf ( "Calculating particle mass via integration with %d grid cells per side\n", N );

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
