#include "thing.h"

template<class T>
std::vector<T*> Thing::getRelation()
{
   std::vector<T*> rels;

   for(auto i: relations) {
      if(T* r0 = dynamic_cast<T*>(i)) {
        rels.push_back(r0);
      }
   }

   return rels;
}

template<class T>
std::vector<T*> Thing::getRelatedObject()
{
   std::vector<T*> objs;

   for(auto i: relations) {
      if(T* r0 = dynamic_cast<T*>(i->getRange())) {
        objs.push_back(r0);
      }
   }

   return objs;
}
