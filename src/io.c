#include "globals.h"
#include "io.h"
#include <inttypes.h>

/* These two handle the F90 Records required for the file format */
#define WRITE_F90REC  {my_fwrite(&blocksize,sizeof(int),1,fp);}

struct GADGET_Header Header;
int blksize;

void Write_output ( bool verbose )
{
    FILE *fp;

    size_t nBytes = 256 + 0.5 * Param.Npart * sizeof ( *P )
                    + 0.5 * Param.Npart * sizeof ( *SphP );


    if ( verbose == 1 ) {
        printf ( "Output : \n"
                 "   File Name = %s\n"
                 "   File Size ~ %.1f MB\n"
                 , Problem.Name, nBytes / 1e6 );
    }

    if ( ! ( fp = fopen ( Problem.Name, "w" ) ) )
        fprintf ( stderr, "Can't open file %s\n"
                  , Problem.Name );

    write_header ( fp, verbose );

    for ( int iblock = 0; iblock < IO_LASTENTRY; iblock++ ) {
        add_block ( fp, ( enum iofields ) iblock, verbose );
    }

    fclose ( fp );

    if ( verbose == 1 ) {
        printf ( "done\n" );
    }

    return;
}

void write_header ( FILE *fp, bool verbose )
{
    if ( verbose == 1 ) {
        printf ( "   Header\n" );
    }

    /* Set Header */
    Header.npart[0] = Param.Npart;
    Header.mass[0] = Problem.Mpart;
    Header.npartTotal[0] = Header.npart[0];

    for ( int i = 1; i < 6; i++ ) {
        Header.npart[i] = Header.mass[i] = Header.npartTotal[i] = 0;
    }

    Header.time = 0;
    Header.redshift = 0;
    Header.flag_sfr = 0;
    Header.flag_feedback = 0;
    Header.flag_cooling = 0;
    Header.num_files = 1;
    Header.BoxSize = Problem.Boxsize[0];
    Header.Omega0 = 0;
    Header.OmegaLambda = 0;
    Header.HubbleParam = 0;
    Header.flag_stellarage = 0;
    Header.flag_metals = 0;

    /* Write F90 Unformatted Record */

    size_t blocksize = sizeof ( int ) + 4 * sizeof ( char );

    WRITE_F90REC                    // <- uses "blocksize"

    my_fwrite ( ( void * ) "HEAD", sizeof ( char ), 4, fp );

    size_t nextblock = sizeof ( Header ) + 2 * sizeof ( int );

    my_fwrite ( &nextblock, sizeof ( int ), 1, fp );

    WRITE_F90REC

    blocksize = sizeof ( Header );

    WRITE_F90REC

    my_fwrite ( &Header, blocksize, 1, fp );

    WRITE_F90REC

    return;
}

void add_block ( FILE *fp, enum iofields iblock, bool verbose )
{
    set_block_info ( ( enum iofields ) iblock );

    if ( verbose == 1 ) {
        printf ( "   Block %d (%s)\n", iblock, Block.Name );
    }

    size_t nData = Block.Ntot * Block.Val_per_element * Block.Bytes_per_element;

    void *write_buffer = Malloc ( nData );

    size_t offset = 0;

    int imax = Block.Npart[0];

    size_t ibuf = 0;

    for ( int type = 0; type < 6; type++ ) {

        for ( int ipart = offset; ipart < imax; ipart++ ) {

            fill_write_buffer ( iblock, write_buffer, ipart, ibuf );

            ibuf += Block.Val_per_element;
        }

        offset += Block.Npart[type];

        if ( type < 5 ) {
            imax += Block.Npart[type + 1];
        }
    }

    size_t blocksize = sizeof ( int ) + 4 * sizeof ( char );

    WRITE_F90REC

    my_fwrite ( &Block.Label, sizeof ( char ), 4, fp );

    size_t nextblock = nData + 2 * sizeof ( int );

    my_fwrite ( &nextblock, sizeof ( int ), 1, fp );

    WRITE_F90REC

    blocksize = nData;

    WRITE_F90REC

    my_fwrite ( write_buffer, blocksize, 1, fp );

    blocksize = nData;

    WRITE_F90REC

    free ( write_buffer );

    return;
}

void fill_write_buffer ( enum iofields blocknr, void *wbuf,
                         size_t ipart, size_t ibuf )
{
    int i = 0;

    switch ( blocknr ) {
    case IO_POS:
        for ( i = 0; i < 3; i++ ) {
            ( ( float * ) wbuf ) [ibuf + i] = P[ipart].Pos[i];
        }
        break;
    case IO_VEL:
        for ( i = 0; i < 3; i++ ) {
            ( ( float * ) wbuf ) [ibuf + i] = P[ipart].Vel[i];
        }
        break;
    case IO_ID:
        ( ( unsigned int * ) wbuf ) [ibuf] = P[ipart].ID;
        break;
    case IO_RHO:
        ( ( float * ) wbuf ) [ibuf] = SphP[ipart].Rho;
        break;
    case IO_RHOMODEL:
        ( ( float * ) wbuf ) [ibuf] = SphP[ipart].Rho_Model;
        break;
    case IO_HSML:
        ( ( float * ) wbuf ) [ibuf] = SphP[ipart].Hsml;
        break;
    case IO_U:
        ( ( float * ) wbuf ) [ibuf] = SphP[ipart].U;
        break;
    case IO_BFLD:
        for ( i = 0; i < 3; i++ ) {
            ( ( float * ) wbuf ) [ibuf + i] = SphP[ipart].Bfld[i];
        }
        break;
    default:
        Assert ( 0, "Block not found %d", blocknr );
        break;
    }

    return;
}
void set_block_info ( enum iofields blocknr )
{
    int i = 0;

    for ( i = Block.Ntot = 0; i < 6; i++ ) {
        Block.Npart[i] = 0;
    }

    strncpy ( Block.Name, "                   ", 22 );

    switch ( blocknr ) {
    case IO_POS:
        strncpy ( Block.Label, "POS ", 4 );       /* Has to be 4 Letters */
        strncpy ( Block.Name, "Coordinates", 13 );
        Block.Npart[0] = Param.Npart;

        Block.Val_per_element = 3;
        Block.Bytes_per_element = sizeof ( P[0].Pos[0] );
        break;
    case IO_VEL:
        strncpy ( Block.Label, "VEL ", 4 );
        strncpy ( Block.Name, "Velocities", 11 );
        Block.Npart[0] = Param.Npart;
        Block.Val_per_element = 3;
        Block.Bytes_per_element = sizeof ( P[0].Vel[0] );
        break;
    case IO_ID:
        strncpy ( Block.Label, "ID  ", 4 );
        strncpy ( Block.Name, "ParticleIDs", 13 );
        Block.Npart[0] = Param.Npart;

        Block.Val_per_element = 1;
        Block.Bytes_per_element = sizeof ( unsigned int );
        break;
    case IO_RHO:
        strncpy ( Block.Label, "RHO ", 4 );
        strncpy ( Block.Name, "Density", 16 );
        Block.Npart[0] = Param.Npart;
        Block.Val_per_element = 1;
        Block.Bytes_per_element = sizeof ( SphP[0].Rho );
        break;
    case IO_RHOMODEL:
        strncpy ( Block.Label, "RHOM", 4 );
        strncpy ( Block.Name, "Model Density", 16 );
        Block.Npart[0] = Param.Npart;
        Block.Val_per_element = 1;
        Block.Bytes_per_element = sizeof ( SphP[0].Rho_Model );
        break;
    case IO_HSML:
        strncpy ( Block.Label, "HSML", 4 );
        strncpy ( Block.Name, "SmoothingLength", 16 );
        Block.Npart[0] = Param.Npart;
        Block.Val_per_element = 1;
        Block.Bytes_per_element = sizeof ( SphP[0].Hsml );
        break;
    case IO_U:
        strncpy ( Block.Label, "U   ", 4 );
        strncpy ( Block.Name, "InternalEnergy", 16 );
        Block.Npart[0] = Param.Npart;
        Block.Val_per_element = 1;
        Block.Bytes_per_element = sizeof ( SphP[0].U );
        break;
    case IO_BFLD:
        strncpy ( Block.Label, "BFLD", 4 );
        strncpy ( Block.Name, "MagneticFIeld", 16 );
        Block.Npart[0] = Param.Npart;
        Block.Val_per_element = 3;
        Block.Bytes_per_element = sizeof ( SphP[0].Bfld[0] );
        break;
    case IO_LASTENTRY:
        strncpy ( Block.Label, "LAST", 4 );
        strncpy ( Block.Name, " ", 1 );
        for ( i = 0; i < 6; i++ ) {
            Block.Npart[i] = 0;
        }
        Block.Bytes_per_element = 0;
        Block.Val_per_element = 0;
        break;
    default:
        fprintf ( stderr, "Block not found %d", blocknr );
        exit ( 5 );
        break;
    }

    for ( int i = Block.Ntot = 0; i < 6; i++ ) {
        Block.Ntot += Block.Npart[i];
    }

    return;
}

/* Basic routine to write data to a file */

size_t my_fwrite ( void *data, size_t size, size_t nItems, FILE *fp )
{
    size_t nWritten = 0;

    if ( nItems * size > 0 ) {

        if ( ( nWritten = fwrite ( data, size, nItems, fp ) ) != nItems ) {

            fprintf ( stderr, "I/O error (fwrite) " );

            fflush ( stderr );

            exit ( 6 );
        }

    } else {
        nWritten = 0;
    }

    return nWritten;
}

#undef WRITE_F90REC

/* Read a number of tags from an ascii file the comment sign is % */

#define LASTPARAMETERID -1111
#define REAL 1
#define STRING 2
#define INT 3

void Read_param_file ( char *filename )
{
    FILE *fd;

    int id[MAXTAGS];
    void *addr[MAXTAGS];
    char tag[MAXTAGS][50];

    char buf[CHARBUFSIZE], buf1[CHARBUFSIZE];
    char buf2[CHARBUFSIZE], buf3[2 * CHARBUFSIZE];
    int tagDone[MAXTAGS] = {0};
    int i, j, nt = 0;

    strcpy ( tag[nt], "Npart" );
    addr[nt] = &Param.Npart;
    id[nt++] = INT;

    strcpy ( tag[nt], "Maxiter" );
    addr[nt] = &Param.Maxiter;
    id[nt++] = INT;

    strcpy ( tag[nt], "MpsFraction" );
    addr[nt] = &Param.MpsFraction;
    id[nt++] = REAL;

    strcpy ( tag[nt], "StepReduction" );
    addr[nt] = &Param.StepReduction;
    id[nt++] = REAL;

    strcpy ( tag[nt], "LimitMps" );
    addr[nt] = &Param.LimitMps[0];
    id[nt++] = REAL;

    strcpy ( tag[nt], "LimitMps10" );
    addr[nt] = &Param.LimitMps[1];
    id[nt++] = REAL;

    strcpy ( tag[nt], "LimitMps100" );
    addr[nt] = &Param.LimitMps[2];
    id[nt++] = REAL;

    strcpy ( tag[nt], "LimitMps1000" );
    addr[nt] = &Param.LimitMps[3];
    id[nt++] = REAL;

    strcpy ( tag[nt], "MoveFraction" );
    addr[nt] = &Param.MoveFraction;
    id[nt++] = REAL;

    strcpy ( tag[nt], "RedistributionFrequency" );
    addr[nt] = &Param.RedistributionFrequency;
    id[nt++] = INT;

    strcpy ( tag[nt], "LastMoveStep" );
    addr[nt] = &Param.LastMoveStep;
    id[nt++] = INT;

    strcpy ( tag[nt], "Problem_Flag" );
    addr[nt] = &Param.Problem_Flag;
    id[nt++] = INT;

    strcpy ( tag[nt], "Problem_Subflag" );
    addr[nt] = &Param.Problem_Subflag;
    id[nt++] = INT;

#ifdef EAT_PNG
    strcpy ( tag[nt], "PNG_Filename" );
    addr[nt] = &Image.Name;
    id[nt++] = STRING;
#endif

    /* Add above */
    id[nt] = LASTPARAMETERID;

    if ( ( fd = fopen ( filename, "r" ) ) ) {
        sprintf ( buf, "%s%s", filename, "-usedvalues" );
        printf ( "\nReading Parameter file : %s \n\n",
                 filename );
        while ( fgets ( buf, CHARBUFSIZE, fd ) ) {
            if ( sscanf ( buf, "%s%s%s", buf1, buf2, buf3 ) < 2 ) {
                continue;
            }

            if ( buf1[0] == '%' ) {
                continue;
            }

            for ( i = 0, j = -1; i < nt; i++ ) {

                if ( ! ( strcmp ( buf1, tag[i] ) ) && ( tagDone[i] != 1 ) ) {
                    j = i;
                    tagDone[i] = 1;
                    break;
                }
            }

            if ( j >= 0 ) {
                switch ( id[j] ) {
                case REAL:
                    * ( ( double * ) addr[j] ) = atof ( buf2 );
                    break;
                case STRING:
                    strcpy ( ( char * ) addr[j], buf2 );
                    break;
                case INT:
                    * ( ( int * ) addr[j] ) = atoi ( buf2 );
                    break;
                }
            }
        }
        fclose ( fd );
    } else {
        fprintf ( stderr, "Parameter file %s not found.\n",
                  filename );
        exit ( 1 );
    }

    for ( i = 0; i < nt; i++ ) {
        if ( !tagDone[i] ) {
            fprintf ( stderr,
                      "Value for tag '%s' missing in parameter file '%s'.\n",
                      tag[i], filename );
            exit ( 1 );
        }
    }

    return;
}
#undef REAL
#undef STRING
#undef INT
#undef LASTPARAMETERID

