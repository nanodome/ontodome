#ifndef CONNECT_H
#define CONNECT_H

#include "baseclass.h"

class Object;

class Connect : public BaseClass
{
    Object* obj0;
    Object* obj1;
public:
    Connect(Object* _obj0, Object* _obj1);

    std::string getClassName() { return "Connect"; }

    Object* getObject0() const;
    Object* getObject1() const;
};

#endif // CONNECT_H
