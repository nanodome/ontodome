#include <iostream>
#include <vector>

#include "base/thing.h"

int main()
{
    Physical universe("universe");

    Item a("H2O");
    Elementary b("H");
    Elementary c("O");
    Elementary d("H");

    ScalarQuantity dens(1,"density","kg/m3");

    addRelation<hasPart,Thing,Thing>(&a,&b);
    addRelation<hasPart,Thing,Thing>(&a,&c);
    addRelation<hasPart,Thing,Thing>(&a,&d);

    addRelation<hasConnection,Thing,Thing>(&b,&c);
    addRelation<hasConnection,Thing,Thing>(&b,&d);
    addRelation<hasConnection,Thing,Thing>(&c,&d);

    addRelation<hasProperty,Thing,Thing>(&a,&dens);

    auto rel = a.getRelation<hasConnection>();

    for(auto &i: rel)
        std::cout << i->getDomain()->getName() << ' '
                  << i->getRelationName() << ' '
                  << i->getRange()->getName() << std::endl;

    return 0;
}
