#include "../globals.h"


float Constant_Density(const int ipart)
{
	double volume = Problem.Boxsize[0]*Problem.Boxsize[1]*Problem.Boxsize[2];
	double mass = Param.Npart * Problem.Mpart;

	return mass/volume;
}
