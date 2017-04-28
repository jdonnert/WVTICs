#include "globals.h"

/* Set IDs with spacing so an ID domain decomposition is more balanced */

void Make_IDs()
{
    printf ( "Make IDs ..." );
    fflush ( stdout );

    #pragma omp parallel for
    for ( int ipart = Param.Npart; ipart < Param.Npart; ipart++ ) {
        P[ipart].ID = ipart + 1;
    }

    size_t delta = 127;

    for ( ;; )
        if ( ( Param.Npart % ++delta ) == 0 ) {
            break;
        }


    printf ( " ID spacing is %zu ...", delta );
    fflush ( stdout );

    int id = 1 - delta, start = 1;

    for ( int ipart = 0; ipart < Param.Npart; ipart++ ) {

        id += delta;

        if ( id > Param.Npart ) {

            start++;

            id = start;
        }

        P[ipart].ID = id;
    }

    printf ( " done\n\n" );
    fflush ( stdout );

    return ;
}
