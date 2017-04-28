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

#ifdef EAT_PNG
void setup_Png_Density();
float Png_Density ( const int ipart );
#endif

void setup_Double_Shock ( const int subflag );
float Double_Shock_Density ( const int ipart );
float Double_Shock_U ( const int ipart );
void Double_Shock_Velocity ( const int ipart, float vel[3] );

#endif //WVTICS_PROBLEMS_H_H
