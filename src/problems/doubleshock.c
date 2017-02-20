#include "../globals.h"

float XBoxhalf;
float Rho[3] = { 0 }, U[3] = { 0 }, Velx[3] = { 0 };

float compression_factor(const float M, const float gamma)
{
	return 1/( 1/p2(M) + (gamma-1)/(gamma+1) * (1 - 1/p2(M)) ); // Mo+ 2010, eq. 8.49
}

float pressure_factor(const float M, const float gamma)
{
	return 2*gamma/(gamma+1)*M*M - (gamma-1)/(gamma+1); // Mo+ 2010, eq. 8.50
}

void Setup_Double_Shock(const int subflag)
{			
	const float gamma = 5.0/3.0;


	
	XBoxhalf = Problem.Boxsize[0]/2;

	float cs = 900; 						// km/s
	float Mach[2] = { 2, 1.5 }; 
  	Rho[0] = 1e-28 * (p3(ULength) / UMass); // rest is Rankine-Hugoniot jump conditions
	
	U[0] = cs*cs/gamma/(gamma-1);
	
	Velx[0] = cs * Mach[0];

	float sigma_v = compression_factor(Mach[0], gamma);
	float sigma_P = pressure_factor(Mach[0], gamma);
	float sigma_T = sigma_P/sigma_v; // Mo+ 2010, eq. 8.51

	Rho[1] = Rho[0] * sigma_v;
	Velx[1] = Velx[0] / sigma_v;
	U[1] = U[0] * sigma_T;

	sigma_v = compression_factor(Mach[1], gamma);
	sigma_P = pressure_factor(Mach[1], gamma);
	sigma_T = sigma_P/sigma_v; // Mo+ 2010, eq. 8.51

	Rho[2] = Rho[1] * sigma_v;
	Velx[2] = Velx[1] / sigma_v;
	U[2] = U[1] * sigma_T;

	Problem.Rho_Max = Rho[2]*1.1;

	printf("Rho=%g %g %g, U=%g %g %g, Vx=%g %g %g \n"
			"Mach=%g %g Sigmav=%g %g \n",
			Rho[0], Rho[1], Rho[2], U[0], U[1], U[2], Velx[0], Velx[1], Velx[2], 
			Mach[0], Mach[1], 
			compression_factor(Mach[0], gamma), compression_factor(Mach[1], gamma));

	return ;
}

float Double_Shock_Density(const int ipart)
{
	const float x = P[ipart].Pos[0];

	float rho = 0;

	if (x < XBoxhalf) 
		rho = Rho[0];
	else if (x < 1.5 * XBoxhalf)
		rho = Rho[1];
	else
		rho = Rho[2];

	return rho;
}

float Double_Shock_U(const int ipart)
{
	const float x = P[ipart].Pos[0];

	float u = 0;

	if (x < XBoxhalf)
		u = U[0];
	else if (x < 1.5 * XBoxhalf)
		u = U[1];
	else
		u = U[2];
	
	
	return u;
}

void Double_Shock_Velocity(const int ipart, float vel[3])
{
	const float x = P[ipart].Pos[0];

	vel[1] = vel[2] = 0;

	if (x < XBoxhalf)
		vel[0] = Velx[0];
	else if (x < 1.5 * XBoxhalf)
		vel[0] = Velx[1];
	else
		vel[0] = Velx[2];

	return ;
}
