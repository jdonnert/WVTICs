#include "globals.h"

double G = 0;

void Set_units()
{
    Unit.Time = Unit.Length/Unit.Vel;
    Unit.Energy = Unit.Mass*Unit.Vel*Unit.Vel;
    Unit.Density = Unit.Mass/p3(Unit.Length);

    printf("Setting System of Units: \n"
            "   Unit Length = %g cm \n"
            "   Unit Time   = %g sec\n"
            "   Unit Mass   = %g g  \n"
            "   Unit Vel    = %g cm/s\n"
            "   Unit Density= %g g/cm^3\n"
            "   Unit Energy = %g erg\n\n"
            ,Unit.Length,Unit.Time, Unit.Mass
            , Unit.Vel, Unit.Density, Unit.Energy);

    G = Grav / p3(Unit.Length) * Unit.Mass * p2(Unit.Time); // internal units

    return;
}

double U2T(double U)
{
    return( (adiabatic_index -1) * U * p2(Unit.Vel) * m_p *
            mean_mol_weight / k_B);
}
double T2U(double T)
{
    return(T/( (adiabatic_index -1) * Unit.Vel * Unit.Vel * m_p *
                mean_mol_weight / k_B));
}

double Density(float rho)
{
    return( ((double)rho) * Unit.Mass / p3(Unit.Length) );
}

