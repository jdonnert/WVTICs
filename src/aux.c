#include "globals.h"

/* Memory Management */

void *Malloc_info ( const char *func, const char *file, const int line,
                    size_t size )
{
    void *result = malloc ( size );

    Assert_Info ( file, func, line, result != NULL || size == 0,
                  "Allocation failed, %zu Bytes \n" , size );

    return result;
}

void *Realloc_info ( const char *func, const char *file, const int line,
                     void *ptr, size_t size )
{
    void *result = realloc ( ptr, size );

    Assert_Info ( file, func, line, result != NULL || size == 0,
                  "Reallocation failed: %zu bytes \n" , size );

    return result;
}

void Free_info ( const char *func, const char *file, const int line, void *ptr )
{
    if ( ptr != NULL ) {
        free ( ptr );
    } else
        fprintf ( stderr, "\nWARNING ! Tried to free a NULL pointer at "
                  "file %s, function %s : line %d \n",
                  file, func, line );
    return ;
}

/* Error Handling, we use variable arguments to be able
 * to print more informative error messages */

void Assert_Info ( const char *func, const char *file, int line,
                   int64_t expr, const char *errmsg, ... )
{
    if ( expr ) {
        return;
    }

    va_list varArgList;

    va_start ( varArgList, errmsg );

    /* we fucked up, tell them */
    fprintf ( stderr,
              "\nERROR : In file %s, function %s(), line %d :\n\n    ",
              file, func, line );

    vfprintf ( stderr, errmsg, varArgList );

    fprintf ( stderr, "\n\n" );

    fflush ( stderr );

    va_end ( varArgList );

    exit ( EXIT_FAILURE ); // ... fatality

    return;
}
