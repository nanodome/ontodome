#ifndef THING_H
#define THING_H

#include <iostream>
#include <vector>
#include <map>

#include "baseclass.h"
#include "relation.h"
#include "datatypes.h"

class Thing : public BaseClass {
//class Thing {

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

    template<class T, class T0>
    void createRelationsTo(std::vector<T0*> o1) {

      for (auto i : o1) {
        T* r = new T(this,i);

        this->addRelation(r);
        i->addRelation(r);
      }
    }

    // find a relation of a specific type
    // TODO: optimize search using multimap<std::string,Relation*>,
    //       getClassName() function for the key and a static_cast<T*>
    template<class T>
    std::vector<T*> getRelations();

    template<class T>
    std::vector<T*> getRelatedObjects();

    template<class T>
    std::vector<double> getRelatedScalarObjects();

    template<class T>
    std::vector<std::vector<double>> getRelatedVectorObjects();
};

class Item : public Thing {
public:
    std::string getClassName() const { return "Item"; }
};

class Collection : public Thing {
public:
    std::string getClassName() const { return "Collection"; }
};

class Physical : public Item {
public:
    std::string getClassName() const { return "Physical"; }
};

class Quantum : public Item {
public:
    std::string getClassName() const { return "Quantum"; }
};

class Void : public Item {
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
    Vector(std::vector<double> s) {data = s;}
    std::string getClassName() const { return "Vector"; }
};

class IUPAC : public String
{
public:
    IUPAC(std::string _iupac) { data = _iupac; }

    std::string getClassName() const { return "IUPAC Name"; }
};

class Unit : public String
{
public:
    Unit(std::string s) {data = s;}
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
    std::string getClassName() const { return "Scalar"; }
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

class PressureTimeDerivative : public ScalarQuantity
{
public:
    PressureTimeDerivative(Scalar* _s, Unit* _u) : ScalarQuantity(_s,_u) {}

    std::string getClassName() const { return "PressureTimeDerivative"; }
};

class Temperature : public ScalarQuantity
{
public:
    Temperature(Scalar* _s, Unit* _u) : ScalarQuantity(_s,_u) {}

    std::string getClassName() const { return "Temperature"; }
};

class TemperatureTimeDerivative : public ScalarQuantity
{
public:
    TemperatureTimeDerivative(Scalar* _s, Unit* _u) : ScalarQuantity(_s,_u) {}

    std::string getClassName() const { return "TemperatureTimeDerivative"; }
};

class MolarFraction : public ScalarQuantity
{
public:
    MolarFraction(Scalar* _s, Unit* _u) : ScalarQuantity(_s,_u) {}

    std::string getClassName() const { return "MolarFraction"; }
};

class Mass : public ScalarQuantity
{
public:
    Mass(Scalar* _s, Unit* _u) : ScalarQuantity(_s,_u) {}

    std::string getClassName() const { return "Mass"; }
};

class BulkDensityLiquid : public ScalarQuantity
{
public:
    BulkDensityLiquid(Scalar* _s, Unit* _u) : ScalarQuantity(_s,_u) {}

    std::string getClassName() const { return "BulkDensityLiquid"; }
};

class BulkDensitySolid : public ScalarQuantity
{
public:
    BulkDensitySolid(Scalar* _s, Unit* _u) : ScalarQuantity(_s,_u) {}

    std::string getClassName() const { return "BulkDensitySolid"; }
};

class MeltingPoint : public ScalarQuantity
{
public:
    MeltingPoint(Scalar* _s, Unit* _u) : ScalarQuantity(_s,_u) {}

    std::string getClassName() const { return "MeltingPoint"; }
};

class Viscosity : public ScalarQuantity
{
public:
    Viscosity(Scalar* _s, Unit* _u) : ScalarQuantity(_s,_u) {}

    std::string getClassName() const { return "Viscosity"; }
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

class SurfaceTension : public ScalarQuantity
{
public:
    SurfaceTension(Scalar* _s, Unit* _u) : ScalarQuantity(_s,_u) {}

    std::string getClassName() const { return "SurfaceTension"; }
};

class SaturationPressure : public ScalarQuantity
{
public:
    SaturationPressure(Scalar* _s, Unit* _u) : ScalarQuantity(_s,_u) {}

    std::string getClassName() const { return "SaturationPressure"; }
};

class KnowledgeGenerator : virtual public Perspective {
public:
    std::string getClassName() const { return "KnowledgeGenerator"; }
};

class Model : virtual public Perspective {
public:
    std::string getClassName() const { return "Model"; }
};

class SoftwareModel : public Model, public KnowledgeGenerator {
public:
    std::string getClassName() const { return "SoftwareModel"; }

    virtual void run() = 0;
};

class MathematicalModel : public Model {
public:
    std::string getClassName() const { return "MathematicalModel"; }
};

class PhysicsBasedModel : public MathematicalModel {
public:
    std::string getClassName() const { return "PhysicsBasedModel"; }
};

class ContinuumModel : public PhysicsBasedModel {
public:
    std::string getClassName() const { return "ContinuumModel"; }
};

class MesoscopicModel : public PhysicsBasedModel {
public:
    std::string getClassName() const { return "MesoscopicModel"; }
};

class Reductionistic : public Perspective {
public:
    std::string getClassName() const { return "Reductionistic"; }
};

class Existent : public Reductionistic {
public:
    std::string getClassName() const { return "Existent"; }
};

class State : public Existent {
public:
    std::string getClassName() const { return "State"; }
};

class Matter : public Perspective {
public:
    std::string getClassName() const { return "Matter"; }
};

class Continuum : public Matter {
public:
    std::string getClassName() const { return "Continuum"; }
};

class Fluid : public Continuum {
public:
    std::string getClassName() const { return "Fluid"; }
};

class Gas : public Fluid {
public:
    std::string getClassName() const { return "Gas"; }
};

class GasMixture : public Gas {
public:
    std::string getClassName() const { return "GasMixture"; }
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
