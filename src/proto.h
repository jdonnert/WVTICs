/* profiles */

double Gas_Density_Profile(const double, const double, const double);
void Read_param_file(char *);
void Set_units();
void Setup();
void Make_positions();
void Make_IDs();
void Make_velocities();
void Make_temperatures();
void Find_sph_quantities();
void Wvt_relax();
void Write_output();
void Write_positions(char*);
void Bfld_from_rotA_SPH();
void Shift_Origin();
void Regularise_sph_particles();
void Show_mass_in_r200();
void Smooth_SPH_quantities();

/* Helper Monkeys */

void *Malloc_info(const char* func, const char* file, const int line,
        size_t size);
void *Realloc_info(const char* func, const char* file, const int line,
        void *ptr, size_t size);
void Free_info(const char* func, const char* file, const int line, void *ptr);
void Assert_Info(const char *func, const char *file, int line, int64_t expr,
        const char *errmsg, ...);

double U2T(double U);
double T2U(double T);
double Density(float rho);

void Print_Int_Bits128(const peanoKey val);
void Print_Int_Bits128r(const peanoKey val);

/* From system libs */
double erand48(unsigned short *);

/* Problem Functions */

float Constant_Density(const int ipart);
float Magneticum_Density(const int ipart);


