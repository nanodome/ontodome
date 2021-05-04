#include "thing.h"

template<class T>
std::vector<T*> Thing::getRelation()
{
  std::vector<T*> rels;
  std::string rel_type;


       for(auto i: relations) {
         if(T* r0 = dynamic_cast<T*>(i)) {
             rels.push_back(r0);
         }
       }
       if (rels.empty()) {
         T* r0 = dynamic_cast<T*>(r0);
         std::string str = typeid(r0).name();
         std::cout << "No " << str.substr(3,str.size()) << " relations found" << std::endl;
         abort();
       }
       else {
         return rels;
       }
}
