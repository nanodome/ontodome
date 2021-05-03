#include <iostream>
#include <vector>

class Thing;

class Relation {
   Thing *o0;
   Thing *o1;
public:
   Relation(Thing* _o0, Thing* _o1)
       : o0(_o0), o1(_o1) {}

   Thing* getDomain() { return o0; }
   Thing* getRange()  { return o1; }

   virtual void getName() = 0;

   virtual ~Relation() = default;
};

class Thing {
   std::vector<Relation*> relations;
public:

   void addRelation(Relation *t)
   {
       relations.push_back(t);
   }

   template<class T>
   auto getRelation()
   {
     std::vector<T*> rels;
     for(auto i: relations)
       if(T* r0 = dynamic_cast<T*>(i)){
           rels.push_back(r0);
         }
     if (rels.empty())
       std::cout << "No relations found" << std::endl;
     return rels;
   }
};

class Object : public Thing {};

class Atom : public Object {};

class hasPart : public Relation {
public:
   hasPart(Object* _o0, Object* _o1) : Relation(_o0,_o1) {}

   void getName() { std::cout << "hasPart" << std::endl; }
};

class hasProperPart : public hasPart {};

class hasConnection : public Relation {};

class hasSign : public Relation {
public:
 hasSign(Object* _o0, Object* _o1) : Relation(_o0,_o1) {}

 void getName() { std::cout << "hasSign" << std::endl; }
};

class hasProperty : public hasSign {
public:
 hasProperty(Object* _o0, Object* _o1) : hasSign(_o0,_o1) {}

 void getName() { std::cout << "hasProperty" << std::endl; }
};

template<class T, class T0, class T1>
Relation* addRelation(T0* o0, T1* o1) {

   T* r = new T(o0,o1);

   o0->addRelation(r);
   o1->addRelation(r);

   return r;
}
