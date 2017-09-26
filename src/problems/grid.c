#include "../globals.h"

float Grid_Density ( const int ipart )
{
    const float x = P[ipart].Pos[0] * Grid.Xpix / Problem.Boxsize[0];
    const float y = P[ipart].Pos[1] * Grid.Ypix / Problem.Boxsize[1];
    const float z = P[ipart].Pos[2] * Grid.Zpix / Problem.Boxsize[1];

    long int index = 1L * floor ( z ) * Grid.Xpix * Grid.Ypix + floor ( y ) * Grid.Xpix + floor ( x );

    return Grid.Density[index];
}
