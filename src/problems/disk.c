#include "../globals.h"

#define Mgas 3.6
#define diskscalelength 2.3
#define diskheigth 0.2

void setup_Exponential_Disk()
{
    Problem.Periodic[0] = Problem.Periodic[1] = Problem.Periodic[2] = 0;

    Problem.Boxsize[0] = 50.0;
    Problem.Boxsize[1] = 50.0;
    Problem.Boxsize[2] = 20.0;

    sprintf ( Problem.Name, "IC_Disk" );

    const double rho = 9 * Mgas / ( 4 * pi * diskscalelength * diskscalelength * diskheigth ); // This value is empiric

    Problem.Rho_Max = rho;

    Density_Func_Ptr = &Exponential_Disk_Density;

}


/* At first we set up a constant density in the Box */
float Exponential_Disk_Density ( const int ipart )
{
    double const x = P[ipart].Pos[0] - Problem.Boxsize[0] * 0.5;
    double const y = P[ipart].Pos[1] - Problem.Boxsize[1] * 0.5;
    double const z = P[ipart].Pos[2] - Problem.Boxsize[2] * 0.5;

    double Radius = sqrt ( x * x + y * y );

    //printf("Mgas: %g\n", Mgas);

    const double cutoffRadius = 10.0, cutoffHeight = 5.0;
    //const double cutoffRadius = 5.0, cutoffHeight = 2.5;
    const double norm = Mgas / ( 4 * pi * diskscalelength * diskscalelength * diskheigth );

    if ( Radius <= cutoffRadius * diskscalelength && fabs ( z ) <= cutoffHeight * diskheigth ) {
        return norm * exp ( -Radius / diskscalelength ) * pow ( 2 / ( exp ( z / diskheigth ) + exp ( -z / diskheigth ) ), 2 );
    } else {
        return norm * exp ( -cutoffRadius ) * pow ( 2 / ( exp ( cutoffHeight ) + exp ( -cutoffHeight ) ), 2 );
    }

}
