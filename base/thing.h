#ifndef THING_H
#define THING_H

#include <iostream>
#include <vector>

#include "relation.h"

class Thing {
    // hardcoded annotations
    std::string name;

    // dynamic relations
    std::vector<Relation*> relations;
public:

    Thing(std::string _name) : name(_name) {}

    virtual std::string getClassName() const { return "Thing"; }

    // annotations
    std::string getName() const { return name; }

    // populate with relations
    void addRelation(Relation *t) { relations.push_back(t); }

    // find a relation of a specific type
    // TODO: optimize search using multimap<std::string,Relation*>,
    //       getClassName() function for the key and a static_cast<T*>
    template<class T>
    std::vector<T*> getRelation();
};

class Item : public Thing {
public:
    Item(std::string _name) : Thing(_name) {}
};

class Collection : public Thing {};

class Physical : public Item {
public:
    Physical(std::string _name) : Item(_name) {}
};

class Quantum : public Item {
public:
    Quantum(std::string _name) : Item(_name) {}
};

class Void : public Item {
public:
    Void(std::string _name) : Item(_name) {}
};

class Elementary : public Physical {
public:
    Elementary(std::string _name) : Physical(_name) {}
};

class Perspective : public Physical {
public:
    Perspective(std::string _name) : Physical(_name) {}
};

#include "thing.cpp"

#endif // THING_H
