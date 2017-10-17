#ifndef IO_H
#define IO_H

struct Block_Info {
    char Name[CHARBUFSIZE];
    char Label[4];
    long long Npart[6];
    long long Ntot;
    long long Ntot_thisHalo;
    size_t Val_per_element;
    size_t Bytes_per_element;
} Block;

// Gadget-2 expects specific block order
enum iofields {
    IO_POS,
    IO_VEL,
    IO_ID,
    IO_U,
    IO_RHO,
    IO_HSML,
    IO_BFLD,
    IO_RHOMODEL,
    IO_LASTENTRY
};

int find_block ( FILE *, char * );

void write_header ( FILE *, bool );
void add_block ( FILE *, enum iofields, bool );
void set_block_info ( enum iofields );
void fill_write_buffer ( enum iofields, void *, size_t, size_t );
size_t  my_fwrite ( void *, size_t, size_t, FILE * );

size_t my_fread ( FILE *file, void *ptr, size_t size, size_t nmemb );
bool checkFileAtEndSafely ( FILE *file );
void swapNbyte ( char *data, int n, int m );

void readGriddedoData();
void writeGridTestData();

void readGadget2File ( const char *fileName );

#endif
