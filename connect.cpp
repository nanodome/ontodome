#include "connect.h"

Connect::Connect(Object* _obj0, Object* _obj1) :
    obj0(_obj0), obj1(_obj1) {}

Object* Connect::getObject0() const { return obj0; }

Object* Connect::getObject1() const { return obj1; }

