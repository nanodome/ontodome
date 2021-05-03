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

    addRelation<hasPart,Thing,Thing>(&a,&b);
    addRelation<hasPart,Thing,Thing>(&a,&c);
    addRelation<hasPart,Thing,Thing>(&a,&d);

    std::vector<hasOverlap*> rel = a.getRelation<hasOverlap>();

    for(auto &i: rel)
        std::cout << i->getDomain()->getName() << ' '
                  << i->getRelationName() << ' '
                  << i->getRange()->getName() << std::endl;
    std::cout << "commit test" << std::endl;

    return 0;
}
