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

class Quantity : public Thing
{
    std::string unit;
public:
    Quantity(std::string _name, std::string _unit) : Thing(_name), unit(_unit) {}

    std::string getClassName() const { return "Quantity"; }

    std::string getUnit() const { return unit; }
};

class ScalarQuantity : public Quantity
{
    double value;
public:
    ScalarQuantity(double _value, std::string _name, std::string _unit) : Quantity(_name, _unit), value(_value) {}

    std::string getClassName() const { return "ScalarQuantity"; }

    double getValue() const { return value; }
};

class VectorQuantity : public Quantity
{
    std::vector<double> value;
public:
    VectorQuantity(const std::vector<double> _value, std::string _name, std::string _unit) : Quantity(_name, _unit), value(_value) {}

    std::string getClassName() const { return "VectorQuantity"; }

    std::vector<double> getValue() const { return value; }
};

class Item : public Thing {
public:
    Item(std::string _name) : Thing(_name) {}

    std::string getClassName() const { return "Item"; }
};

class Collection : public Thing {
public:
    Collection(std::string _name) : Thing(_name) {}

    std::string getClassName() const { return "Collection"; }
};

class Physical : public Item {
public:
    Physical(std::string _name) : Item(_name) {}

    std::string getClassName() const { return "Physical"; }
};

class Quantum : public Item {
public:
    Quantum(std::string _name) : Item(_name) {}

    std::string getClassName() const { return "Quantum"; }
};

class Void : public Item {
public:
    Void(std::string _name) : Item(_name) {}

    std::string getClassName() const { return "Void"; }
};

class Elementary : public Physical {
public:
    Elementary(std::string _name) : Physical(_name) {}

    std::string getClassName() const { return "Elementary"; }
};

class Perspective : public Physical {
public:
    Perspective(std::string _name) : Physical(_name) {}

    std::string getClassName() const { return "Perspective"; }
};

#include "thing.cpp"

#endif // THING_H
