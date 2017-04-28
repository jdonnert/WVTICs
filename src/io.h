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
#ifdef SPH_CUBIC_SPLINE
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
#else
enum iofields {
    IO_POS,
    IO_VEL,
    IO_ID,
    IO_RHO,
    IO_RHOMODEL,
    IO_HSML,
    IO_U,
    IO_BFLD,
    IO_LASTENTRY
};
#endif  // SPH_CUBIC_SPLINE

size_t my_fread ( void *, size_t, size_t, FILE *, int );
int find_block ( FILE *, char * );

void write_header ( FILE *, bool );
void add_block ( FILE *, enum iofields, bool );
void set_block_info ( enum iofields );
void fill_write_buffer ( enum iofields, void *, size_t, size_t );
size_t  my_fwrite ( void *, size_t, size_t, FILE * );
