#include "thing.h"

template<class T>
std::vector<T*> Thing::getRelations()
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
std::vector<T*> Thing::getRelatedObjects()
{
   std::vector<T*> objs;

   for(auto i: relations) {
      if(T* r0 = dynamic_cast<T*>(i->getRange())) {
        objs.push_back(r0);
      }
   }

   return objs;
}

template<class T>
std::vector<double> Thing::getRelatedScalarObjects()
{
   std::vector<double> objs;

   for(auto i: relations) {
      if(T* r0 = dynamic_cast<T*>(i->getRange())) {
        objs.push_back(r0->template getRelatedObjects<Scalar>()[0]->data);
      }
   }

   return objs;
}

template<class T>
std::vector<std::vector<double>> Thing::getRelatedVectorObjects()
{
   std::vector<std::vector<double>> objs;

   for(auto i: relations) {
      if(T* r0 = dynamic_cast<T*>(i->getRange())) {
        objs.push_back(r0->template getRelatedObjects<Vector>()[0]->data);
      }
   }

   return objs;
}
