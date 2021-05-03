#include "thing.h"

template<class T>
std::vector<T*> Thing::getRelation()
{
    std::vector<T*> out;

    for(auto i: relations)
        if(T* r0 = dynamic_cast<T*>(i))
             out.push_back(r0);

    return out;
}
