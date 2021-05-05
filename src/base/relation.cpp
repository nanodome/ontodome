#include "relation.h"

template<class T, class T0, class T1>
Relation* addRelation(T0* o0, T1* o1) {

    T* r = new T(o0,o1);

    o0->addRelation(r);
    o1->addRelation(r);

    return r;
}

template<class T, class T0, class T1>
std::vector<Relation*> addRelations(T0* o0, std::vector<T1*> o1) {

  std::vector<Relation*> rr;
  for (auto i : o1){
    T* r = new T(o0,i);

    o0->addRelation(r);
    i->addRelation(r);
  }
    return rr;
}
