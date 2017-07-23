#include "../globals.h"
#include "ryu_jones.h"

struct Ryu_Jones_Parameters RJ95 = { 0 };

void setup_Ryu_Jones_Shocktube()
{
    printf ( "Ryu Jones Shocktube\n"
             "   parameter\tleft \t\tright \n"
             "   density \t%f \t%f \n"
             "   vx      \t%f \t%f \n"
             "   vy      \t%f \t%f \n"
             "   vz      \t%f \t%f \n"
             "   bx      \t%f \t%f \n"
             "   by      \t%f \t%f \n"
             "   bz      \t%f \t%f \n"
             "   energy  \t%f \t%f \n\n",
             RJ95.left_rho, RJ95.right_rho, RJ95.left_vx, RJ95.right_vx,
             RJ95.left_vy,  RJ95.right_vy, RJ95.left_vz, RJ95.right_vz,
             RJ95.left_bx, RJ95.right_bx, RJ95.left_by, RJ95.right_by,
             RJ95.left_bz, RJ95.right_bz, RJ95.left_energy, RJ95.right_energy );

    Problem.Boxsize[0] = 140.0;
    Problem.Boxsize[1] = 10.0;
    Problem.Boxsize[2] = 10.0;

    const double rho_mean = 0.5 * ( RJ95.left_rho + RJ95.right_rho );

    Problem.Rho_Max = max ( RJ95.left_rho, RJ95.right_rho );
    Problem.Mpart = rho_mean * ( Problem.Boxsize[0] * Problem.Boxsize[1] * Problem.Boxsize[2] ) / Param.Npart;

    Density_Func_Ptr = &Ryu_Jones_Density;
    U_Func_Ptr = &Ryu_Jones_U;
    Velocity_Func_Ptr = &Ryu_Jones_Velocity;
    Magnetic_Field_Func_Ptr = &Ryu_Jones_Magnetic_Field;
    return;
}

bool isOnLeft ( const int ipart )
{
    const double x = P[ipart].Pos[0];

    if ( x <= 0.5 * Problem.Boxsize[0] ) {
        return true;
    } else {
        return false;
    }
}

float Ryu_Jones_Density ( const int ipart )
{
    if ( isOnLeft ( ipart ) ) {
        return RJ95.left_rho;
    } else {
        return RJ95.right_rho;
    }
}

void Ryu_Jones_Velocity ( const int ipart, float out[3] )
{
    if ( isOnLeft ( ipart ) ) {
        out[0] = RJ95.left_vx;
        out[1] = RJ95.left_vy;
        out[2] = RJ95.left_vz;
    } else {
        out[0] = RJ95.right_vx;
        out[1] = RJ95.right_vy;
        out[2] = RJ95.right_vz;
    }
}

void Ryu_Jones_Magnetic_Field ( const int ipart, float out[3] )
{
    if ( isOnLeft ( ipart ) ) {
        out[0] = RJ95.left_bx;
        out[1] = RJ95.left_by;
        out[2] = RJ95.left_bz;
    } else {
        out[0] = RJ95.right_bx;
        out[1] = RJ95.right_by;
        out[2] = RJ95.right_bz;
    }
}

float Ryu_Jones_U ( const int ipart )
{
    if ( isOnLeft ( ipart ) ) {
        return RJ95.left_energy;
    } else {
        return RJ95.right_energy;
    }
}

void set_Ryu_Jones_Parameters_1A()
{
    sprintf ( Problem.Name, "IC_Ryu_Jones_1A" );

    // e.g. https://www.astro.umn.edu/groups/compastro/?q=node/1521
    RJ95.left_rho = 1;
    RJ95.left_vx = 10;
    RJ95.left_vy = 0;
    RJ95.left_vz = 0;
    RJ95.left_bx = 5 / sqrt ( 4 * pi );
    RJ95.left_by = 5 / sqrt ( 4 * pi );
    RJ95.left_bz = 0;
    RJ95.left_energy = 20;

    RJ95.right_rho = 1;
    RJ95.right_vx = -10;
    RJ95.right_vy = 0;
    RJ95.right_vz = 0;
    RJ95.right_bx = 5 / sqrt ( 4 * pi );
    RJ95.right_by = 5 / sqrt ( 4 * pi );
    RJ95.right_bz = 0;
    RJ95.right_energy = 1;
}

void set_Ryu_Jones_Parameters_1B()
{
    sprintf ( Problem.Name, "IC_Ryu_Jones_1B" );

    // e.g. https://www.astro.umn.edu/groups/compastro/?q=node/1531
    RJ95.left_rho = 1;
    RJ95.left_vx = 0;
    RJ95.left_vy = 0;
    RJ95.left_vz = 0;
    RJ95.left_bx = 3 / sqrt ( 4 * pi );
    RJ95.left_by = 5 / sqrt ( 4 * pi );
    RJ95.left_bz = 0;
    RJ95.left_energy = 1;

    RJ95.right_rho = 0.1;
    RJ95.right_vx = 0;  // WOMBAT -10
    RJ95.right_vy = 0;
    RJ95.right_vz = 0;
    RJ95.right_bx = 3 / sqrt ( 4 * pi );
    RJ95.right_by = 2 / sqrt ( 4 * pi );
    RJ95.right_bz = 0;
    RJ95.right_energy = 10;
}

void set_Ryu_Jones_Parameters_2A()
{
    sprintf ( Problem.Name, "IC_Ryu_Jones_2A" );

    // e.g. https://www.astro.umn.edu/groups/compastro/?q=node/1541
    RJ95.left_rho = 1.08;
    RJ95.left_vx = 1.2;
    RJ95.left_vy = 0.01;
    RJ95.left_vz = 0.5;
    RJ95.left_bx = 2 / sqrt ( 4 * pi );
    RJ95.left_by = 3.6 / sqrt ( 4 * pi );
    RJ95.left_bz = 2 / sqrt ( 4 * pi );
    RJ95.left_energy = 0.95;

    RJ95.right_rho = 1;
    RJ95.right_vx = 0;
    RJ95.right_vy = 0;
    RJ95.right_vz = 0;
    RJ95.right_bx = 2 / sqrt ( 4 * pi );
    RJ95.right_by = 4 / sqrt ( 4 * pi );
    RJ95.right_bz = 2 / sqrt ( 4 * pi );  // WOMBAT 5 / sqrt ( 4 * pi );
    RJ95.right_energy = 1;
}

void set_Ryu_Jones_Parameters_2B()
{
    sprintf ( Problem.Name, "IC_Ryu_Jones_2B" );

    // e.g. https://www.astro.umn.edu/groups/compastro/?q=node/1551
    RJ95.left_rho = 1;
    RJ95.left_vx = 0;
    RJ95.left_vy = 0;
    RJ95.left_vz = 0;
    RJ95.left_bx = 3 / sqrt ( 4 * pi );
    RJ95.left_by = 6 / sqrt ( 4 * pi );
    RJ95.left_bz = 0;
    RJ95.left_energy = 1;

    RJ95.right_rho = 0.1;
    RJ95.right_vx = 0;
    RJ95.right_vy = 2;  // WOMBAT 0;
    RJ95.right_vz = 1;  // WOMBAT 0;
    RJ95.right_bx = 3 / sqrt ( 4 * pi );
    RJ95.right_by = 1 / sqrt ( 4 * pi );
    RJ95.right_bz = 0;
    RJ95.right_energy = 10;
}

void set_Ryu_Jones_Parameters_3A()
{
    sprintf ( Problem.Name, "IC_Ryu_Jones_3A" );

    // e.g. https://www.astro.umn.edu/groups/compastro/?q=node/1561
    RJ95.left_rho = 0.1;
    RJ95.left_vx = 50;
    RJ95.left_vy = 0;
    RJ95.left_vz = 0;
    RJ95.left_bx = 0;
    RJ95.left_by = -1 / sqrt ( 4 * pi );
    RJ95.left_bz = -2 / sqrt ( 4 * pi );
    RJ95.left_energy = 0.4;

    RJ95.right_rho = 0.1;
    RJ95.right_vx = 0;
    RJ95.right_vy = 0;
    RJ95.right_vz = 0;
    RJ95.right_bx = 0;
    RJ95.right_by = 1 / sqrt ( 4 * pi );
    RJ95.right_bz = 2 / sqrt ( 4 * pi );
    RJ95.right_energy = 0.2;
}

void set_Ryu_Jones_Parameters_3B()
{
    sprintf ( Problem.Name, "IC_Ryu_Jones_3B" );

    // e.g. https://www.astro.umn.edu/groups/compastro/?q=node/1571
    RJ95.left_rho = 1;
    RJ95.left_vx = -1;
    RJ95.left_vy = 0;
    RJ95.left_vz = 0;
    RJ95.left_bx = 0;
    RJ95.left_by = 1;
    RJ95.left_bz = 0;
    RJ95.left_energy = 1;

    RJ95.right_rho = 1;
    RJ95.right_vx = 1;
    RJ95.right_vy = 0;
    RJ95.right_vz = 0;
    RJ95.right_bx = 0;
    RJ95.right_by = 1;
    RJ95.right_bz = 0;
    RJ95.right_energy = 1;
}

void set_Ryu_Jones_Parameters_4A()
{
    sprintf ( Problem.Name, "IC_Ryu_Jones_4A" );

    // e.g. https://www.astro.umn.edu/groups/compastro/?q=node/1581
    RJ95.left_rho = 1;
    RJ95.left_vx = 0;
    RJ95.left_vy = 0;
    RJ95.left_vz = 0;
    RJ95.left_bx = 1;
    RJ95.left_by = 1;
    RJ95.left_bz = 0;
    RJ95.left_energy = 1;

    RJ95.right_rho = 0.2;
    RJ95.right_vx = 0;
    RJ95.right_vy = 0;
    RJ95.right_vz = 0;
    RJ95.right_bx = 1;
    RJ95.right_by = 0;
    RJ95.right_bz = 0;
    RJ95.right_energy = 0.1;
}

void set_Ryu_Jones_Parameters_4B()
{
    sprintf ( Problem.Name, "IC_Ryu_Jones_4B" );

    // e.g. https://www.astro.umn.edu/groups/compastro/?q=node/1591
    RJ95.left_rho = 0.4;
    RJ95.left_vx =  -0.66991;
    RJ95.left_vy = 0.98263;
    RJ95.left_vz = 0;
    RJ95.left_bx = 1.3;
    RJ95.left_by = 0.0025293;
    RJ95.left_bz = 0;
    RJ95.left_energy = 0.52467;

    RJ95.right_rho = 1;
    RJ95.right_vx = 0;
    RJ95.right_vy = 0;
    RJ95.right_vz = 0;
    RJ95.right_bx = 1.3;
    RJ95.right_by = 1;
    RJ95.right_bz = 0;
    RJ95.right_energy = 1;
}

void set_Ryu_Jones_Parameters_4C()
{
    sprintf ( Problem.Name, "IC_Ryu_Jones_4C" );

    // e.g. https://www.astro.umn.edu/groups/compastro/?q=node/1601
    RJ95.left_rho = 0.65;
    RJ95.left_vx = 0.667; // WOMBAT -0.667;
    RJ95.left_vy = -0.257;
    RJ95.left_vz = 0;
    RJ95.left_bx = 0.75;
    RJ95.left_by = 0.55;
    RJ95.left_bz = 0;
    RJ95.left_energy = 0.5;

    RJ95.right_rho = 1;
    RJ95.right_vx = 0.4;
    RJ95.right_vy = -0.94;
    RJ95.right_vz = 0;
    RJ95.right_bx = 0.75;
    RJ95.right_by = 0;
    RJ95.right_bz = 0;
    RJ95.right_energy = 0.75;
}

void set_Ryu_Jones_Parameters_4D()
{
    sprintf ( Problem.Name, "IC_Ryu_Jones_4D" );

    // e.g. https://www.astro.umn.edu/groups/compastro/?q=node/1611
    RJ95.left_rho = 1;
    RJ95.left_vx = 0;
    RJ95.left_vy = 0;
    RJ95.left_vz = 0;
    RJ95.left_bx = 0.7;
    RJ95.left_by = 0;
    RJ95.left_bz = 0;
    RJ95.left_energy = 1;

    RJ95.right_rho = 0.3;
    RJ95.right_vx = 0;
    RJ95.right_vy = 0;
    RJ95.right_vz = 1;
    RJ95.right_bx = 0.7;
    RJ95.right_by = 1;
    RJ95.right_bz = 0;
    RJ95.right_energy = 0.2;
}

void set_Ryu_Jones_Parameters_5A()
{
    sprintf ( Problem.Name, "IC_Ryu_Jones_5A" );

    // e.g. https://www.astro.umn.edu/groups/compastro/?q=node/1621
    RJ95.left_rho = 1;
    RJ95.left_vx = 0;
    RJ95.left_vy = 0;
    RJ95.left_vz = 0;
    RJ95.left_bx = 0.75;
    RJ95.left_by = 1;
    RJ95.left_bz = 0;
    RJ95.left_energy = 1;

    RJ95.right_rho = 0.125;
    RJ95.right_vx = 0;
    RJ95.right_vy = 0;
    RJ95.right_vz = 0;
    RJ95.right_bx = 0.75;
    RJ95.right_by = -1;
    RJ95.right_bz = 0;
    RJ95.right_energy = 0.1;
}

void set_Ryu_Jones_Parameters_5B()
{
    sprintf ( Problem.Name, "IC_Ryu_Jones_5B" );

    // e.g. https://www.astro.umn.edu/groups/compastro/?q=node/1631
    RJ95.left_rho = 1;
    RJ95.left_vx = 0;
    RJ95.left_vy = 0;
    RJ95.left_vz = 0;
    RJ95.left_bx = 1.3;
    RJ95.left_by = 1;
    RJ95.left_bz = 0;
    RJ95.left_energy = 1;

    RJ95.right_rho = 0.4;
    RJ95.right_vx = 0;
    RJ95.right_vy = 0;
    RJ95.right_vz = 0;
    RJ95.right_bx = 1.3;
    RJ95.right_by = -1;
    RJ95.right_bz = 0;
    RJ95.right_energy = 0.4;
}
