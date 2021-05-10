#ifndef RELATION_H
#define RELATION_H

#include <string>
#include <vector>

#include "baseclass.h"

class Thing;

class Relation : public BaseClass {
    Thing *o0;
    Thing *o1;
public:
    Relation(Thing* _o0, Thing* _o1)
        : o0(_o0), o1(_o1) {}

    virtual std::string getRelationName() const { return "Relation"; }

    Thing* getDomain() { return o0; }
    Thing* getRange()  { return o1; }

    virtual ~Relation() = default;
};

class isModelFor : public Relation {
public:
    isModelFor(Thing* _o0, Thing* _o1) : Relation(_o0,_o1) {}
};

class hasOverlap : public Relation {
public:
    hasOverlap(Thing* _o0, Thing* _o1) : Relation(_o0,_o1) {}
};

class hasPart : public hasOverlap {
public:
    hasPart(Thing* _o0, Thing* _o1) : hasOverlap(_o0,_o1) {}

    virtual std::string getRelationName() const { return "hasPart"; }
};

class hasProperPart : public hasPart {
public:
    hasProperPart(Thing* _o0, Thing* _o1) : hasPart(_o0,_o1) {
        if(_o0==_o1)
            abort();
    }

    std::string getRelationName() const { return "hasProperPart"; }
};

class hasConnection : public Relation {
public:
  hasConnection(Thing* _o0, Thing* _o1) : Relation(_o0,_o1) {}

  std::string getRelationName() const { return "hasConnection"; }
};

class hasSign : public Relation {
public:
  hasSign(Thing* _o0, Thing* _o1) : Relation(_o0,_o1) {}
};

class hasProperty : public hasSign {
public:
  hasProperty(Thing* _o0, Thing* _o1) : hasSign(_o0,_o1) {}

  std::string getRelationName() const { return "hasProperty"; }
};

class hasScalarProperty : public hasProperty {
public:
  hasScalarProperty(Thing* _o0, Thing* _o1) : hasProperty(_o0,_o1) {}

  std::string getRelationName() const { return "hasScalarProperty"; }
};



// Relation based functions
template<class T, class T0, class T1>
Relation* addRelation(T0* o0, T1* o1);

// Add multiple relation at once for lazy people
template<class T, class T0, class T1>
std::vector<Relation*> addRelations(T0* o0, std::vector<T1*> o1);

#include "relation.cpp"

#endif // RELATION_H
