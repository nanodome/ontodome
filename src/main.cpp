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

    addRelation<hasProperty,Thing,Thing>(&a,&dens);

    std::vector<hasProperty*> rel = a.getRelation<hasProperty>();

    for(auto &i: rel)
        std::cout << i->getDomain()->getName() << ' '
                  << i->getRelationName() << ' '
                  << i->getRange()->getName() << std::endl;

    return 0;
}
