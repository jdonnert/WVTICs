#include "gadgetreader.h"
#include "../globals.h"

void readHeader ( FILE *file )
{
    int blksize, blocksize, dummysize, buffer;

    blocksize = findBlock ( file, "HEAD" );
    if ( blocksize <= 0 ) {
        fprintf ( stderr, "Block <HEAD> not found!\n" );
        exit ( 5 );
    }

    dummysize = blocksize - 16 * sizeof ( int ) - 12 * sizeof ( double );

    my_fread ( file,  &blksize, sizeof ( int ), 1 );
    //swapNbyte ( ( char* ) &blksize, 1, sizeof ( int ) );

    my_fread ( file,  Header.npart, 6 * sizeof ( int ), 1 );
    //swapNbyte ( ( char* ) Header.npart, 6, sizeof ( int ) );

    my_fread ( file,  Header.mass, 6 * sizeof ( double ), 1 );
    //swapNbyte ( ( char* ) Header.mass, 6, sizeof ( double ) );

    my_fread ( file,  &Header.time, sizeof ( double ), 1 );
    //swapNbyte ( ( char* ) &Header.time, 1, sizeof ( double ) );

    my_fread ( file,  &Header.redshift, sizeof ( double ), 1 );
    //swapNbyte ( ( char* ) &Header.redshift, 1, sizeof ( double ) );

    my_fread ( file,  &buffer, sizeof ( int ), 1 );
    //swapNbyte ( ( char* ) buffer, 1, sizeof ( int ) );
    Header.flag_sfr = buffer;

    my_fread ( file,  &buffer, sizeof ( int ), 1 );
    //swapNbyte ( ( char* ) buffer, 1, sizeof ( int ) );
    Header.flag_feedback = buffer;

    fseek ( file, 6 * sizeof ( int ), 1 );

    my_fread ( file,  &buffer, sizeof ( int ), 1 );
    //swapNbyte ( ( char* ) buffer, 1, sizeof ( int ) );
    Header.flag_cooling = buffer;

    my_fread ( file,  &buffer, sizeof ( int ), 1 );
    //swapNbyte ( ( char* ) buffer, 1, sizeof ( int ) );
    Header.num_files = buffer;

    my_fread ( file,  &Header.BoxSize, sizeof ( double ), 1 );
    //swapNbyte ( ( char* ) &Header.BoxSize, 1, sizeof ( double ) );

    my_fread ( file,  &Header.Omega0, sizeof ( double ), 1 );
    //swapNbyte ( ( char* ) &Header.Omega0, 1, sizeof ( double ) );

    my_fread ( file,  &Header.OmegaLambda, sizeof ( double ), 1 );
    //swapNbyte ( ( char* ) &Header.OmegaLambda, 1, sizeof ( double ) );

    my_fread ( file,  &Header.HubbleParam, sizeof ( double ), 1 );
    //swapNbyte ( ( char* ) &Header.HubbleParam, 1, sizeof ( double ) );

    fseek ( file, dummysize, 1 );

    my_fread ( file,  &blksize, sizeof ( int ), 1 );
    //swapNbyte ( ( char* ) &blksize, 1, 4 );
}

DataPair readArray ( const char *label, const ParticleTypeInfo &info, const Header *const header )
{
    std::cout << "Reading " << label << std::endl << std::flush;
    int blksize, blocksize;

    blocksize = findBlock ( label );
    if ( blocksize <= 0 ) {
        std::cerr << "Block <" << label << "> not found!" << std::endl << std::flush;
        return DataPair();
    } else {
        const InfoBlock::Info *infopart = infos.findBlockInfo ( label );
        int dimension = infopart->dimension;
        if ( header != 0 ) {
            dimension = blocksize / header->ntot ( infopart->types ) / sizeof ( float );
            if ( dimension != infopart->dimension ) {
                std::cerr << "Info supplied wrong dimension for this block (" << label << ") ! Got " << infopart->dimension << " but expected " << dimension << "!" << std::endl << std::flush;
                std::cerr << "Blocksize: " << blocksize << "; npart: " << header->ntot ( infopart->types ) << std::endl << std::flush;
            }
        }
        //! \todo Additionally use types info to check if selected types are valid and how much seek is needed when skipping blocks

#ifdef MY_DEBUG
        std::cout << "Reading " << info.length *dimension *sizeof ( float ) << " bytes of data (into a " << dimension
                  << " dimensional array) from <" << label << ">" << std::endl;
#endif
        float *data = new float[dimension * info.length];

        my_fread ( file,  &blksize, sizeof ( int ), 1 );
        //swapNbyte ( ( char* ) &blksize, 1, 4 );

        size_t last = 0;
        for ( std::vector<ParticleTypeInfo::Doublet>::const_iterator i = info.ranges.begin(); i < info.ranges.end(); i++ ) {
            if ( i->start - last > 0 ) {
                fseek ( file, i->start - last, 1 );
            }

            size_t l = i->end - i->start;
            my_fread ( file,  & ( data[last] ), dimension * sizeof ( float ), l );

            last = i->end;
        }

        //swapNbyte ( ( char* ) data, dimension * info.length, sizeof ( float ) );

        my_fread ( file,  &blksize, sizeof ( int ), 1 );
        //swapNbyte ( ( char* ) &blksize, 1, 4 );

        DataPair ret ( dimension, dimension * info.length );
        unsigned int remaining = dimension * info.length;
        unsigned int offset = 0u;
        for ( unsigned int dc = 0; dc < ret.getDataCount(); ++dc ) {
            GETDATA ( ret, dc );

            unsigned int l = ret.getLength ( dc );
            if ( l > remaining ) {
                l = remaining;
            }
            std::copy ( data + offset, data + offset + l, retLoc );
            remaining -= l;
            offset += l;
        }

        delete[] data;

        return ret;
    }
}

int findBlock ( FILE *file, const char *label )
{
    int blksize, blocksize = 0;
    char blocklabel[5];
    blocklabel[4] = '\0';

    rewind ( file );

    while ( !checkFileAtEndSafely() && blocksize == 0 ) {
        my_fread ( file,  &blksize, sizeof ( int ), 1 );
        //swapNbyte((char *) &blksize, 1, 4);
        if ( blksize == 134217728 ) { //8 with wrong endian
            //swap = 1 - swap;
            //swapNbyte ( ( char* ) &blksize, 1, 4 );
        }
        if ( blksize != 8 ) {
            fprintf ( stdout, "Incorrect format (blksize = %d)!\n", blksize );
            exit ( 1 );
        } else {
            my_fread ( file,  blocklabel, 4 * sizeof ( char ), 1 );
            my_fread ( file,  &blocksize, sizeof ( int ), 1 );
            //swapNbyte((char *) &blocksize, 1, 4);
            my_fread ( file,  &blksize, sizeof ( int ), 1 );
            //swapNbyte((char *) &blksize, 1, 4);

            if ( strcmp ( label, blocklabel ) != 0 ) {
                if ( strcmp ( blocklabel, infos.getLastBlockName() ) == 0 ) { //If last block encountered
                    return -1;
                } else {
                    fseek ( file, blocksize, 1 );
                    blocksize = 0;
                }
            }
        }
    }
    return ( blocksize - 8 );
}
