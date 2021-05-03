#include <iostream>

#include "base/thing.h"

int main()
{
   Object a;
   Atom b;
   Atom c;

   addRelation<hasPart,Object,Object>(&a,&b);

//    a.getRelation<hasPart>()->getName();

   addRelation<hasProperty,Object,Object>(&b,&c);
   addRelation<hasPart,Object,Object>(&c,&b);
   addRelation<hasPart,Object,Object>(&c,&a);

   auto rels = c.getRelation<hasSign>();
//    for (auto i: c.getRelation<hasPart>())
//      std::cout << "Relation type is: " << i->getName() << std::endl;

//    std::cout << "Hello World!" << std::endl;
   return 0;
}
