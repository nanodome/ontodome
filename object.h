#ifndef OBJECT_H
#define OBJECT_H

#include <vector>
#include <map>

#include "baseclass.h"
#include "quantity.h"
#include "connect.h"

class Object : public BaseClass
{
    std::vector<Object*> parts;
    std::vector<Connect*> connections;
    std::map<std::string,Quantity*> properties;
public:
    Object();

    std::string getClassName() { return "Object"; }

    void add(Connect*) {}
    void add(Object*) {}
    void add(Quantity*) {}
};

#endif // OBJECT_H
