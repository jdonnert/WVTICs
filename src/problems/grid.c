#include "../globals.h"
#include "../io.h"

void setup_Grid_Densiy()
{
    sprintf ( Problem.Name, "IC_Grid" );

    readGriddedoData();

    Density_Func_Ptr = &Grid_Density;

    const long int nElements = Grid.Xpix * Grid.Ypix * Grid.Zpix;
    Problem.Rho_Max = Grid.Density[0];
    for ( long int i = 1; i < nElements; ++i ) {
        if ( Problem.Rho_Max < Grid.Density[0] ) {
            Problem.Rho_Max = Grid.Density[0];
        }
    }

    Problem.Boxsize[0] = Grid.Xpix;
    Problem.Boxsize[1] = Grid.Ypix;
    Problem.Boxsize[2] = Grid.Zpix;
}

float Grid_Density ( const int ipart )
{
    const float x = P[ipart].Pos[0] * Grid.Xpix / Problem.Boxsize[0];
    const float y = P[ipart].Pos[1] * Grid.Ypix / Problem.Boxsize[1];
    const float z = P[ipart].Pos[2] * Grid.Zpix / Problem.Boxsize[1];

    long int index = 1L * floor ( z ) * Grid.Xpix * Grid.Ypix + floor ( y ) * Grid.Xpix + floor ( x );

    return Grid.Density[index];
}
