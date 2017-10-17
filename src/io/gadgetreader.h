#ifndef GADGETREADER_H
#define GADGETREADER_H

#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include "../io.h"

extern struct GADGET_Header {
    int npart[6];
    double mass[6];
    double time;
    double redshift;
    int flag_sfr;
    int flag_feedback;
    unsigned int npartTotal[6];
    int flag_cooling;
    int num_files;
    double BoxSize;
    double Omega0;
    double OmegaLambda;
    double HubbleParam;
    int flag_stellarage;
    int flag_metals;
    unsigned int npartTotalHighWord[6];
    char fill[64];      // fill to 256 Bytes
} Header;

void readHeader ( FILE *file );
void readBlock ( FILE *file, enum iofields block );

DataPair readArray ( const char *label, const ParticleTypeInfo &info, const Header *const header = 0 );

int findBlock ( FILE *file, const char *label );

#endif
