#ifndef BASECLASS_H
#define BASECLASS_H

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

boost::uuids::random_generator uuid_gen;

class BaseClass
{
    // uuid
    boost::uuids::uuid uuid;

public:
    BaseClass() { uuid = uuid_gen(); }

    boost::uuids::uuid getUuid() const { return uuid; }

};

#endif // BASECLASS_H
