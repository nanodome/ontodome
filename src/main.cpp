#include <iostream>
#include <vector>

#include "base/thing.h"
#include "nanodome.h"
#include "gasphase/gasphasecv.h"
//#include "species/species.h"

int main()
{
    Physical universe("universe");

    HeteronuclearMolecule a("H2O");
    Atom b("H");
    Atom c("O");
    Atom d("H");

    ScalarQuantity dens(1,"density","kg/m3");
    VectorQuantity pos({1,1,1},"posi","m");

    addRelation<hasPart,Thing,Thing>(&a,&b);
    addRelation<hasPart,Thing,Thing>(&a,&c);
    addRelation<hasPart,Thing,Thing>(&a,&d);

    addRelation<hasProperty,Thing,Thing>(&a,&dens);
    addRelation<hasProperty,Thing,Thing>(&a,&pos);


    auto res = a.getScalarProperty("density");
    auto resp = a.getVectorProperty("posi");

    auto rel = a.getRelation<hasPart>();

    for(auto &i: rel) {
        std::cout << i->getDomain()->getName() << ' '
                  << i->getRelationName() << ' '
                  << i->getRange()->getName() << std::endl;
    }

    HeteronuclearMolecule si("Si");

    // Create the Species Object with the given properties
    ScalarQuantity Mass(28.085*AMU,"mass","kg");
    ScalarQuantity T_Melt(1687,"melting point","K");
    ScalarQuantity Sigma(3.30e-10,"L-J sigma","m");
    ScalarQuantity Eps(4.37e-20,"L-J epsilon","J");
    ScalarQuantity Bulk_Dens_Liq(2570,"liquid bulk density","kg/m3");
    ScalarQuantity Bulk_Dens_Sol(2329,"solid bulk density","kg/m3");
    VectorQuantity S_Ten({0.732, 0.000086, 1685.},"surface tension coefficients","#");
    VectorQuantity P_Sat({7.5341, 23399.},"saturation pressure coefficients","#");
    ScalarQuantity Visc(5e-4,"viscosity","Pa*s");

    addRelations<hasProperty,Thing,Thing>(&si,{&Mass,&T_Melt,&Sigma,&Eps,&Bulk_Dens_Liq,&Bulk_Dens_Sol,&S_Ten,&P_Sat,&Visc});

    HeteronuclearMolecule fe("Fe");

    // Create the Species Object with the given properties
    ScalarQuantity Mass1(140*AMU,"mass","kg");
    ScalarQuantity T_Melt1(1687,"melting point","K");
    ScalarQuantity Sigma1(3.30e-10,"L-J sigma","m");
    ScalarQuantity Eps1(4.37e-20,"L-J epsilon","J");
    ScalarQuantity Bulk_Dens_Liq1(2570,"liquid bulk density","kg/m3");
    ScalarQuantity Bulk_Dens_Sol1(2329,"solid bulk density","kg/m3");
    VectorQuantity S_Ten1({0.732, 0.000086, 1685.},"surface tension coefficients","#");
    VectorQuantity P_Sat1({7.5341, 23399.},"saturation pressure coefficients","#");
    ScalarQuantity Visc1(5e-5,"viscosity","Pa*s");

    addRelations<hasProperty,Thing,Thing>(&fe,{&Mass1,&T_Melt1,&Sigma1,&Eps1,&Bulk_Dens_Liq1,&Bulk_Dens_Sol1,&S_Ten1,&P_Sat1,&Visc1});

    std::vector<MolecularEntity*> species = {&si,&fe};
    std::valarray<double> cc = {0.5,0.5};

    GasPhaseCV gp(101325,300,species,cc);
    std::cout << gp.get_viscosity() << std::endl;
    std::cout << gp.get_average_molecular_mass() << std::endl;

    return 0;
}
