#ifndef OBJECT_H
#define OBJECT_H

#include <vector>
#include <map>
#include <iostream>
#include <boost/algorithm/cxx11/any_of.hpp>

#include "baseclass.h"
#include "quantity.h"
#include "connect.h"

class Object : public BaseClass
{
    std::vector<Object*> parts;
//    std::vector<Connect*> connections;
    std::vector<Object*> connections;
    std::map<std::string,ScalarQuantity*> propertiesS;
    std::map<std::string,VectorQuantity*> propertiesV;
    std::string name; //could be IUPAC formula or anything else

public:
    Object(std::string _name);

    std::string getClassName() { return "Object"; }

    /// Methods that represent hasPart and hasProperty relations
//    void add(Connect *_connection);
    void connect(Object *_obj);
    void add(Object *_obj);
    void add(ScalarQuantity *_scalar);
    void add(VectorQuantity *_vector);

    /// Method to get a specified property
    double getScalarProperty(std::string name){
      if (propertiesS.find(name) != propertiesS.end()){
        return propertiesS.at(name)->getValue();
        }
      else {
          std::cout << "The requested property " << name << " is not available. Please add it." << std::endl;
          abort();
        }
    }

    std::valarray<double> getVectorProperty(std::string name){
      if (propertiesV.find(name) != propertiesV.end()){
        return propertiesV.at(name)->getValue();
        }
      else {
          std::cout << "The requested property " << name << " is not available. Please add it." << std::endl;
          abort();
        }
    }
};

#endif // OBJECT_H
