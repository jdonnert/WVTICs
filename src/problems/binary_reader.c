#include "../globals.h"

float *densityBuffer;
long filelen, dimlen;

void setup_BinaryReader_Density ( const int subflag )
{
    Problem.Boxsize[0] = 1.0;
    Problem.Boxsize[1] = 1.0;
    Problem.Boxsize[2] = 1.0;

    Problem.Periodic[0] = true;

    sprintf ( Problem.Name, "IC_Binary" );

    static FILE *infile;
    int error = 0;

    // @todo use own parameter for filename, do not rely on PNG structure
    if ( ! ( infile = fopen ( Image.Name, "rb" ) ) ) {
        fprintf ( stderr, PROG_NAME ":  can't open PNG file [%s]\n", Image.Name );
        ++error;
        exit ( 1 );
    } else {
        fseek ( infile, 0, SEEK_END );
        filelen = ftell ( infile );
        rewind ( infile );

        densityBuffer = ( float * ) malloc ( filelen * sizeof ( char ) );
        fread ( densityBuffer, filelen, 1, infile );
    }
    fclose ( infile );

    dimlen = round ( pow ( filelen / sizeof ( char ), 1.0 / 3.0 ) );

    Problem.Rho_Max = densityBuffer[0];
    for ( long i = 0; i < filelen / sizeof ( char ); ++i ) {
        Problem.Rho_Max = fmax ( Problem.Rho_Max, densityBuffer[i] );
    }

    Density_Func_Ptr = &BinaryReader_Density;
}

float BinaryReader_Density ( const int ipart )
{
    const int x = P[ipart].Pos[0] / Problem.Boxsize[0] * dimlen;
    const int y = P[ipart].Pos[1] / Problem.Boxsize[1] * dimlen;
    const int z = P[ipart].Pos[2] / Problem.Boxsize[2] * dimlen;

    int index = x * dimlen * dimlen + y * dimlen + z;

    return densityBuffer[index];
}
