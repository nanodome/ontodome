#include <iostream>
#include <vector>

#include "base/thing.h"
#include "nanodome.h"
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
    ScalarQuantity T_Melt(1687,"Melting point","K");
    ScalarQuantity Sigma(3.30e-10,"L-J sigma","m");
    ScalarQuantity Eps(4.37e-20,"L-J epsilon","J");
    ScalarQuantity Bulk_Dens_Liq(2570,"Liquid Bulk Density","kg/m3");
    ScalarQuantity Bulk_Dens_Sol(2329,"Solid Bulk Density","kg/m3");
    VectorQuantity S_Ten({0.732, 0.000086, 1685.},"Surface tensione coefficients","#");
    VectorQuantity P_Sat({7.5341, 23399.},"Saturation Pressure Coefficients","#");

    addRelation<hasProperty,Thing,Thing>(&si, &Mass);
    addRelation<hasProperty,Thing,Thing>(&si, &T_Melt);
    addRelation<hasProperty,Thing,Thing>(&si, &Sigma);
    addRelation<hasProperty,Thing,Thing>(&si, &Eps);
    addRelation<hasProperty,Thing,Thing>(&si, &Bulk_Dens_Liq);
    addRelation<hasProperty,Thing,Thing>(&si, &Bulk_Dens_Sol);
    addRelation<hasProperty,Thing,Thing>(&si, &S_Ten);
    addRelation<hasProperty,Thing,Thing>(&si,&P_Sat);

    return 0;
}
