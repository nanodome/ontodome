#ifndef RELATION_H
#define RELATION_H

#include <string>

class Thing;

class Relation {
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
};

class hasConnection : public Relation {};

class hasSign : public Relation {};

class hasProperty : public hasSign {};


// Relation based functions
template<class T, class T0, class T1>
Relation* addRelation(T0* o0, T1* o1);

#include "relation.cpp"

#endif // RELATION_H
