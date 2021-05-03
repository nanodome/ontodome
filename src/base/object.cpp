#include "object.h"

Object::Object(std::string _name) : name(_name) {}

//void Object::add(Connect *_connection)
//{
//  connections.push_back(_connection);
//  std::cout << "Added connection: " << _connection << std::endl;
//}

void Object::connect(Object *_obj)
{
  connections.push_back(_obj);

  std::cout << "Added connection to: " << _obj << std::endl;
  if (!boost::algorithm::any_of_equal(_obj->connections, this)){
    _obj->connect(this);
    }
}

void Object::add(Object *_obj)
{
  parts.push_back(_obj);
  std::cout << "Added part: " << _obj << std::endl;
}

void Object::add(ScalarQuantity *_scalar)
{
  propertiesS.insert({_scalar->getName(), _scalar});
  std::cout << "Added scalar quantity: " << _scalar << std::endl;
}

void Object::add(VectorQuantity *_vector)
{
  propertiesV.insert({_vector->getName(), _vector});
  std::cout << "Added vector quantity: " << &_vector << std::endl;
}
