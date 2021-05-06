#include <iostream>
#include <vector>

#include "base/thing.h"

int main()
{
    Pressure p(new Scalar(2), new Unit("atm"));

    Physical obj;

    obj.createRelationTo<hasProperty,Pressure>(&p);

    std::cout << obj.getRelatedObject<Pressure>()[0]->getRelatedObject<Scalar>()[0]->data << std::endl;
    std::cout << obj.getRelatedObject<Pressure>()[0]->getRelatedObject<Unit>()[0]->data << std::endl;
}
