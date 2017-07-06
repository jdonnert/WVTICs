#ifndef WVTICS_PROBLEMS_H_H
#define WVTICS_PROBLEMS_H_H

void setup_Constant_Density();
float Constant_Density ( const int ipart );

void setup_TopHat_Density();
float TopHat_Density ( const int ipart );

void setup_Sawtooth_Density();
float Sawtooth_Density ( const int ipart );

void setup_SineWave_Density();
float SineWave_Density ( const int ipart );

void setup_Gradient_Density();
float Gradient_Density ( const int ipart );

void setup_Magneticum_Density();
float Magneticum_Density ( const int ipart );

void setup_Png_Density();
float Png_Density ( const int ipart );

void setup_Double_Shock ( const int subflag );
float Double_Shock_Density ( const int ipart );
float Double_Shock_U ( const int ipart );
void Double_Shock_Velocity ( const int ipart, float vel[3] );

void setup_Sod_Shock();
float Sod_Shock_Density ( const int ipart );
float Sod_Shock_U ( const int ipart );

void setup_Sedov_Blast();
float Sedov_Blast_Density ( const int ipart );
float Sedov_Blast_U ( const int ipart );
void Sedov_Blast_PostProcessing();

void setup_Kelvin_Helmholtz_Instability();
float Kelvin_Helmholtz_Instability_Density ( const int ipart );
float Kelvin_Helmholtz_Instability_U ( const int ipart );
void Kelvin_Helmholtz_Instability_Velocity ( const int ipart, float out[3] );

void setup_Gresho_Vortex();
float Gresho_Vortex_Density ( const int ipart );
void Gresho_Vortex_Velocity ( const int ipart, float out[3] );
float Gresho_Vortex_U ( const int ipart );

void setup_Blob();
float Blob_Density ( const int ipart );
void Blob_Velocity ( const int ipart, float out[3] );
float Blob_U ( const int ipart );

void setup_Zeldovich_Pancake();
float Zeldovich_Pancake_Density ( const int ipart );
void Zeldovich_Pancake_Velocity ( const int ipart, float out[3] );
float Zeldovich_Pancake_U ( const int ipart );

void setup_Box();
bool isInnerBox ( const int ipart );
float Box_Density ( const int ipart );
void Box_Velocity ( const int ipart, float out[3] );
float Box_U ( const int ipar );

void setup_Rotor();
float Rotor_Density ( const int ipart );
void Rotor_Velocity ( const int ipart, float out[3] );
void Rotor_Magnetic_Field ( const int ipart, float out[3] );
float Rotor_U ( const int ipart );

void setup_Evrard_Collapse();
float Evrard_Collapse_Density ( const int ipart );
void Evrard_Collapse_Velocity ( const int ipart, float out[3] );
float Evrard_Collapse_U ( const int ipart );

void setup_Orszag_Tang_Vortex();
float Orszag_Tang_Vortex_Density ( const int ipart );
void Orszag_Tang_Vortex_Velocity ( const int ipart, float out[3] );
void Orszag_Tang_Vortex_Magnetic_Field ( const int ipart, float out[3] );
float Orszag_Tang_Vortex_U ( const int ipart );

void setup_User_Density ( const int subflag );
float User_Density ( const int ipart );
float User_U ( const int ipart );
void User_Velocity ( const int ipart, float out[3] );
void User_Bfld ( const int ipart, float out[3] );
void User_PostProcessing();

#endif //WVTICS_PROBLEMS_H_H
