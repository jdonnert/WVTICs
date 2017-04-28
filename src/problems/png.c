#ifdef EAT_PNG

#include "../globals.h"
#include "../readpng.h"

void setup_Png_Density()
{
#ifndef EAT_PNG
    printf ( "Error: must use OPT += -DEAT_PNG in Makefile\n" );
    exit ( 1 );
#else

    sprintf ( Problem.Name, "IC_PNG" );

    Setup_Density_From_Image();

    Density_Func_Ptr = &Png_Density;

    Problem.Boxsize[0] = Image.Xpix;
    Problem.Boxsize[1] = Image.Ypix;
    Problem.Boxsize[2] = Image.Zpix;
#endif // EAT_PNG
}

/*
 * Read x,y density from 2D grayscale png image
 */

void Setup_Density_From_Image()
{
    /* http://www.libpng.org/pub/png/book/chapter13.html */
    readpng_version_info();

    static FILE *infile;
    static unsigned long xpix, ypix, image_rowbytes;
    int rc;
    int error = 0;

    if ( ! ( infile = fopen ( Image.Name, "rb" ) ) ) {
        fprintf ( stderr, PROG_NAME ":  can't open PNG file [%s]\n", Image.Name );
        ++error;
        exit ( 1 );
    } else {
        if ( ( rc = readpng_init ( infile, &xpix, &ypix ) ) != 0 ) {
            switch ( rc ) {
            case 1:
                fprintf ( stderr, PROG_NAME
                          ":  [%s] is not a PNG file: incorrect signature\n",
                          Image.Name );
                break;
            case 2:
                fprintf ( stderr, PROG_NAME
                          ":  [%s] has bad IHDR (libpng longjmp)\n", Image.Name );
                break;
            case 4:
                fprintf ( stderr, PROG_NAME ":  insufficient memory\n" );
                break;
            default:
                fprintf ( stderr, PROG_NAME
                          ":  unknown readpng_init() error\n" );
                break;
            }
            ++error;
        }
        if ( error ) {
            fclose ( infile );
        }
    }

    printf ( "\nxpix, ypix = %lu, %lu\n", xpix, ypix );

    static unsigned char bg_red, bg_green, bg_blue = 0;

    if ( readpng_get_bgcolor ( &bg_red, &bg_green, &bg_blue ) > 1 ) {
        readpng_cleanup ( TRUE );
        fprintf ( stderr, PROG_NAME
                  ":  libpng error while checking for background color\n" );
        exit ( 2 );
    }

    printf ( "\nbg r,g,b = %d, %d, %d\n", bg_red, bg_green, bg_blue );

    static int image_channels;
    static unsigned char *image_data;

    image_data = readpng_get_image ( 0, &image_channels, &image_rowbytes );

    if ( !image_data ) {
        fprintf ( stderr, PROG_NAME ":  unable to decode PNG image\n" );
        exit ( 3 );
    }

    printf ( "Beginning data loop (image_channels == %d)\n", image_channels );
    printf ( "(width = %ld, rowbytes = %ld, image_rowbytes = %ld)\n",
             xpix, ypix, image_rowbytes );

    int red, green, blue;

    // Allocate a flat aray to hold the density
    Image.Xpix = ( int ) xpix;
    Image.Ypix = ( int ) ypix;
    Image.Zpix = ( int ) ypix;
    int nBytes = Image.Xpix * Image.Ypix *  sizeof ( float );
    printf ( "xpix, ypix, nBytes = %d %d %d\n", Image.Xpix, Image.Ypix, nBytes );
    Image.Density = Malloc ( nBytes );
    memset ( Image.Density, 0, nBytes );

    float gray;
    unsigned char *src;
    float rho;
    int index;
    for ( int row = 0;  row < Image.Ypix;  row++ ) {
        src = image_data + row * image_rowbytes;
        if ( image_channels == 3 ) {
            for ( int i = 0; i < Image.Xpix;  i++ ) {
                red   = *src++;
                green = *src++;
                blue  = *src++;
                // YUV color modell: rgb -> luminosity
                gray = ( 0.2989 * red + 0.5870 * green + 0.1140 * blue );
                rho = ( 255 - gray ) / 255.;
                // printf("gray = %f; rho = %f\n", gray, rho);
                index = row * Image.Xpix + i;
                Image.Density[index] = rho;
                // printf("Image.Density[%d] = %g for gray = %g\n",
                //    index, Image.Density[index], gray);
            }
        }
    }

    printf ( "Success building Image.Density\n" );

    return;
}

float Png_Density ( const int ipart )
{
    const float x = P[ipart].Pos[0] * Image.Xpix / Problem.Boxsize[0];
    const float y = P[ipart].Pos[1] * Image.Ypix / Problem.Boxsize[1];

    int index = ( floor ( Image.Ypix - y ) * Image.Xpix + floor ( x ) );

    // printf("x, y, z, index, density = %g %g %g, %d, %g\n", x, y, z, index,  Image.Density[index]);
    return Image.Density[index];
}

#endif  // EAT_PNG
