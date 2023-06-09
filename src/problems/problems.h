#ifndef WVTICS_PROBLEMS_H_H
#define WVTICS_PROBLEMS_H_H

void setup_Constant_Density();
float Constant_Density ( const int ipart , const double bias );

void setup_TopHat_Density();
float TopHat_Density ( const int ipart , const double bias );

void setup_Sawtooth_Density();
float Sawtooth_Density ( const int ipart , const double bias );

void setup_SineWave_Density();
float SineWave_Density ( const int ipart , const double bias );

void setup_Gradient_Density();
float Gradient_Density ( const int ipart , const double bias );

void setup_Magneticum_Density();
float Magneticum_Density ( const int ipart , const double bias );

void setup_Png_Density();
float Png_Density ( const int ipart , const double bias);

void setup_Double_Shock ( const int subflag );
float Double_Shock_Density ( const int ipart , const double bias );
float Double_Shock_U ( const int ipart );
void Double_Shock_Velocity ( const int ipart, float vel[3] );

void setup_Sod_Shock();
float Sod_Shock_Density ( const int ipart , const double bias );
float Sod_Shock_U ( const int ipart );

void setup_Sedov_Blast();
float Sedov_Blast_Density ( const int ipart , const double bias);
float Sedov_Blast_U ( const int ipart );
float Sedov_Blast_abs ();
float Sedov_Blast_kinetic ();
void Sedov_Blast_PostProcessing();

void setup_Kelvin_Helmholtz_Instability();
float Kelvin_Helmholtz_Instability_Density ( const int ipart , const double bias );
float Kelvin_Helmholtz_Instability_U ( const int ipart );
void Kelvin_Helmholtz_Instability_Velocity ( const int ipart, float out[3] );

void setup_Keplerian_Ring();
float Keplerian_Ring_Density ( const int ipart , const double bias);
float Keplerian_Ring_U ( const int ipart );
void Keplerian_Ring_Velocity ( const int ipart, float out[3] );

void setup_Blob();
float Blob_Density ( const int ipart , const double bias );
void Blob_Velocity ( const int ipart, float out[3] );
float Blob_U ( const int ipart );

void setup_Evrard_Collapse();
float Evrard_Collapse_Density ( const int ipart , const double bias );
void Evrard_Collapse_Velocity ( const int ipart, float out[3] );
float Evrard_Collapse_U ( const int ipart );

void setup_Zeldovich_Pancake();
float q_of_x ( const int ipart );
float Zeldovich_Pancake_Density ( const int ipart , const double bias );
void Zeldovich_Pancake_Velocity ( const int ipart, float out[3] );
float Zeldovich_Pancake_U ( const int ipart );

void setup_Box();
bool isInnerBox ( const int ipart );
float Box_Density ( const int ipart , const double bias );
void Box_Velocity ( const int ipart, float out[3] );
float Box_U ( const int ipar );

void setup_Gresho_Vortex();
float Gresho_Vortex_Density ( const int ipart , const double bias );
void Gresho_Vortex_Velocity ( const int ipart, float out[3] );
float Gresho_Vortex_U ( const int ipart );

void setup_Rotor();
float Rotor_Density ( const int ipart , const double bias );
void Rotor_Velocity ( const int ipart, float out[3] );
void Rotor_Magnetic_Field ( const int ipart, float out[3] );
float Rotor_U ( const int ipart );

void setup_Strong_Blast();
float Strong_Blast_Density ( const int ipart , const double bias);
void Strong_Blast_Velocity ( const int ipart, float out[3] );
void Strong_Blast_Magnetic_Field ( const int ipart, float out[3] );
float Strong_Blast_U ( const int ipart );

void setup_Exponential_Disk();
float Exponential_Disk_Density ( const int ipart , const double bias );

void setup_Boss();
float Boss_Phi ( double const x, double const y );
float Boss_Density ( const int ipart , const double bias);

void setup_Sphere();
float Sphere_Density ( const int ipart , const double bias );

void setup_Orszag_Tang_Vortex();
float Orszag_Tang_Vortex_Density ( const int ipart , const double bias);
void Orszag_Tang_Vortex_Velocity ( const int ipart, float out[3] );
void Orszag_Tang_Vortex_Magnetic_Field ( const int ipart, float out[3] );
float Orszag_Tang_Vortex_U ( const int ipart );

void setup_Linear_Alfven_Wave();
float Linear_Alfven_Wave_Density ( const int ipart , const double bias );
void Linear_Alfven_Wave_Velocity ( const int ipart, float out[3] );
void Linear_Alfven_Wave_Magnetic_Field ( const int ipart, float out[3] );
float Linear_Alfven_Wave_U ( const int ipart );

void setup_Rayleigh_Taylor_Instability();
float Rayleigh_Taylor_Instability_Density ( const int ipart , const double bias );
void Rayleigh_Taylor_Instability_Velocity ( const int ipart, float out[3] );
void Rayleigh_Taylor_Instability_Magnetic_Field ( const int ipart, float out[3] );
float Rayleigh_Taylor_Instability_U ( const int ipart );

void set_Ryu_Jones_Parameters_1A();
void set_Ryu_Jones_Parameters_1B();
void set_Ryu_Jones_Parameters_2A();
void set_Ryu_Jones_Parameters_2B();
void set_Ryu_Jones_Parameters_3A();
void set_Ryu_Jones_Parameters_3B();
void set_Ryu_Jones_Parameters_4A();
void set_Ryu_Jones_Parameters_4B();
void set_Ryu_Jones_Parameters_4C();
void set_Ryu_Jones_Parameters_4D();
void set_Ryu_Jones_Parameters_5A();
void set_Ryu_Jones_Parameters_5B();
void setup_Ryu_Jones_Shocktube();
bool isOnLeft ( const int ipart );
float Ryu_Jones_Density ( const int ipart, const double bias );
void Ryu_Jones_Velocity ( const int ipart, float out[3] );
void Ryu_Jones_Magnetic_Field ( const int ipart, float out[3] );
float Ryu_Jones_U ( const int ipart );

void setup_User_Density ( const int subflag );
float User_Density ( const int ipart , const double bias );
float User_U ( const int ipart );
void User_Velocity ( const int ipart, float out[3] );
void User_Bfld ( const int ipart, float out[3] );
void User_PostProcessing();

void setup_GalaxyCluster();

#endif //WVTICS_PROBLEMS_H_H
