#include "../globals.h"

static const float gamma = 5.0/3.0; // adiabatic index
static float XBoxhalf;
static float Rho[3] = { 0 }, U[3] = { 0 }, Velx[3] = { 0 };

static double Params[][3] = {
							 {850, 2, 1.6 * 2000},
						 	 {850, 3, 2.4 * 2000},
							 {850, 4, 2.9 * 2000}
};

float compression_factor(const float M, const float gamma)
{
	return 1/( 1/p2(M) + (gamma-1)/(gamma+1) * (1 - 1/p2(M)) ); // Mo+ 2010, eq. 8.49
}

float pressure_factor(const float M, const float gamma)
{
	return 2*gamma/(gamma+1)*M*M - (gamma-1)/(gamma+1); // Mo+ 2010, eq. 8.50
}


float find_M1(const float cs_up, const float v_dw, const float gamma)
{
	double left = 1;
	double right = 100;
	double M = 0;

	for (;;) { // bisection (Press+ 1992)
	
		M = left + 0.5 * (right-left);

		double sigma = compression_factor(M, gamma);
		double res =  sigma*M - v_dw/cs_up; // residual
		
		if (fabs(res) < 1e-4)
			break;

		if (res < 0)
			left = M;
		else
			right = M;
	}

	return M;
}

void set_shock_tube(float cs[3], float Mach[3]) {
	
	Velx[0] = 0; // frame of upstream gas

	/* First Shock */

	float sigma_v = compression_factor(Mach[0], gamma); // Mo+ 2010, eq. 8.51
	float sigma_P = pressure_factor(Mach[0], gamma);
	float sigma_T = sigma_P/sigma_v; 

	Rho[1] = Rho[0] * sigma_v;
	Velx[1] = cs[0] * Mach[0] * (1 - 1/sigma_v);
	U[1] = U[0] * sigma_T;
	cs[1] = sqrt(U[1] * gamma * (gamma-1));

	/* Second Shock */
	Mach[1] = find_M1(cs[1], Velx[2], gamma);

	sigma_v = compression_factor(Mach[1], gamma);
	sigma_P = pressure_factor(Mach[1], gamma);
	sigma_T = sigma_P/sigma_v;

	Rho[2] = Rho[1] * sigma_v;
	U[2] = U[1] * sigma_T;
	cs[2] = sqrt(U[2]*gamma * (gamma-1));

	return;
}

void Setup_Double_Shock(const int subflag)
{			

	XBoxhalf = Problem.Boxsize[0]/2;

	float cs[3] = { Params[subflag][0] ,0 ,0 }; // km/s
	float Mach[2] = { Params[subflag][1], 0 }; 
	Velx[2] = Params[subflag][2];

	Rho[0] = 1e-28 * (p3(ULength) / UMass);
	U[0] = cs[0]*cs[0]/gamma/(gamma-1);

	set_shock_tube(cs, Mach);
	
	Problem.Rho_Max = Rho[2]*1.1;

	printf("Flag %d\n   Rho = (%g %g %g) \n   U = (%g %g %g) \n   Vx = (%g %g %g) \n"
			"   cs = (%g %g %g) \n   Mach = (%g %g) \n   Sigmav = (%g %g) \n\n", 
			subflag, Rho[0], Rho[1], Rho[2], U[0], U[1], U[2], Velx[0], Velx[1], Velx[2], 
			cs[0], cs[1], cs[2],Mach[0], Mach[1], 
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
