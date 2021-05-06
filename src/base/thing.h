#ifndef THING_H
#define THING_H

#include <iostream>
#include <vector>
#include <boost/variant/variant.hpp>
#include <variant>

#include "baseclass.h"
#include "relation.h"
#include "datatypes.h"

class Thing : BaseClass {

    // dynamic relations
    std::vector<Relation*> relations;
public:

    virtual std::string getClassName() const { return "Thing"; }

    // populate with relations
    void addRelation(Relation *t) { relations.push_back(t); }

    template<class T, class T0>
    void createRelationTo(T0* o1) {

        T* r = new T(this,o1);

        this->addRelation(r);
        o1->addRelation(r);
    }

    // find a relation of a specific type
    // TODO: optimize search using multimap<std::string,Relation*>,
    //       getClassName() function for the key and a static_cast<T*>
    template<class T>
    std::vector<T*> getRelation();

    template<class T>
    std::vector<T*> getRelatedObject();
};


class Item : public Thing {
public:
    std::string getClassName() const { return "Item"; }
};

class Collection : public Thing {
public:
    std::string getClassName() const { return "Collection"; }
};

class Physical : public Thing {
public:
    std::string getClassName() const { return "Physical"; }
};

class Quantum : public Item {
public:
    std::string getClassName() const { return "Quantum"; }
};

class Void : public Thing {
public:
    std::string getClassName() const { return "Void"; }
};

class Elementary : public Physical {
public:
    std::string getClassName() const { return "Elementary"; }
};

class Perspective : public Physical {
public:
    std::string getClassName() const { return "Perspective"; }
};

class Symbolic : public Perspective {
public:
    std::string getClassName() const { return "Symbolic"; }
};

class String : public Symbolic, public DataType<std::string>
{
public:
    std::string getClassName() const { return "String"; }
};

class Vector : public Symbolic, public DataType<std::vector<double>>
{
public:
    std::string getClassName() const { return "Vector"; }
};

class Unit : public String
{
public:
    Unit(std::string _unit) { data = _unit; }

    std::string getClassName() const { return "Unit"; }
};

class Quantity : public Symbolic
{
public:
    std::string getClassName() const { return "Quantity"; }
};

class Scalar : public Symbolic, public DataType<double>
{
public:
    Scalar(double s) {data = s;}
    std::string getClassName() const { return "Quantity"; }
};


class ScalarQuantity : public Quantity
{
public:
    ScalarQuantity(Scalar* _s, Unit* _u)
    {
        createRelationTo<hasPart,Thing>(_u);
        createRelationTo<hasPart,Thing>(_s);
    }

    std::string getClassName() const { return "ScalarQuantity"; }
};

class Pressure : public ScalarQuantity
{
public:
    Pressure(Scalar* _s, Unit* _u) : ScalarQuantity(_s,_u) {}

    std::string getClassName() const { return "Pressure"; }
};

class VectorQuantity : public Quantity
{
public:
    VectorQuantity(Vector* _s, Unit* _u)
    {
        createRelationTo<hasPart,Thing>(_u);
        createRelationTo<hasPart,Thing>(_s);
    }

    std::string getClassName() const { return "VectorQuantity"; }
};

class Model : public Perspective {
public:
    std::string getClassName() const { return "Model"; }
};

class PhysicsBasedModel : public Model {
public:
    std::string getClassName() const { return "PhysicsBasedModel"; }
};

class ContinuumModel : public PhysicsBasedModel {
public:
    std::string getClassName() const { return "ContinuumModel"; }
};

class Matter : public Perspective {
public:
    std::string getClassName() const { return "Matter"; }
};

class MolecularEntity : public Matter {
public:
    std::string getClassName() const { return "MolecularEntity"; }
};

class Atom : public MolecularEntity {
public:
    std::string getClassName() const { return "Atom"; }
};

class PolyatomicEntity : public MolecularEntity {
public:
    std::string getClassName() const { return "PolyatomicEntity"; }
};

class HeteronuclearMolecule : public PolyatomicEntity {
public:
    std::string getClassName() const { return "HeteronuclearMolecule"; }
};

class HomonuclearMolecule : public PolyatomicEntity {
public:
    std::string getClassName() const { return "HomonuclearMolecule"; }
};

#include "thing.cpp"

#endif // THING_H
