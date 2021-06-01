#include "thing.h"

template<class T>
std::vector<T*> Thing::getRelations() {
   std::vector<T*> rels;

   for(auto i: relations) {
      if(T* r0 = dynamic_cast<T*>(i)) {
        rels.push_back(r0);
      }
   }

   return rels;
}

// avoids looping on same object by returning the other object in a given relation
Thing* avoider (Thing* o0, Relation* o1) {
  Thing* tested;
  if (o1->getRange()->getUuid() != o0->getUuid()) {
      tested = o1->getRange();
  } else {
      tested = o1->getDomain();
  }
  return tested;
}

// add only UUID to list if not already present
void adder(std::vector<boost::uuids::uuid>* scans, boost::uuids::uuid id) {
  if (std::find(scans->begin(), scans->end(), id) == scans->end()) {
    scans->push_back(id);
  }
}

// looks into an objects relations set and decides wheter to continue looking in its related objects or not
template<class T>
T* Thing::looker(Thing* start, std::vector<boost::uuids::uuid>* scanned) {
  T* result = nullptr;
  adder(scanned,start->getUuid());
  for (auto i : start->relations) {
    auto tested = avoider(start,i);
    if (std::find(scanned->begin(), scanned->end(), tested->getUuid()) == scanned->end()) {
      adder(scanned,tested->getUuid()); // add to scanned list
      if (dynamic_cast<T*>(tested)) {
        result = dynamic_cast<T*>(tested);
        return result;
        break;
      }
    }
  }

  if (result == nullptr) {
    for (auto i : start->relations) {
      auto tested = avoider(start,i);
      adder(scanned,tested->getUuid()); // add to scanned list

      // look in tested object's relations, if there is an unscanned object, mark tested to be looked.
      bool to_scan = false;
      for (auto ii : tested->relations) {
          if (std::find(scanned->begin(), scanned->end(), avoider(tested,ii)->getUuid()) == scanned->end()) {
              to_scan = true;
              break;
          }
      }

      if (to_scan) {
        result = looker<T>(tested, scanned);
        if (result != nullptr) { break; }
      }
    }
  }
  return result;
}

template<class T>
T* Thing::find() {
  auto* start = this; // the starting object will always be the objects from which the method is called
  T* result = nullptr;
  std::vector<boost::uuids::uuid> scanned; // list of scanned objects by means of UUID
  adder(&scanned,start->getUuid()); // add start to scanned items in order to avoid looping

  // scan for object of type T in start object relations
  for (auto i : start->relations) {
    auto tested = avoider(start,i);
    adder(&scanned,tested->getUuid()); //add to scanned list
    if (dynamic_cast<T*>(tested)) {
      result = dynamic_cast<T*>(tested);
      break; // stop searching, we got it
    }
  }

  // if not found above, try searching in related objects
  if (result == nullptr) {
    for (auto i : start->relations) {
      auto tested = avoider(start,i);
      adder(&scanned,tested->getUuid()); // add to scanned list
      result = looker<T>(tested, &scanned);
      if (result != nullptr) { break; } // stop searching, we got it
    }
  }

  if (result == nullptr) {
      std::cout << "Not found." << std::endl;
      abort(); }
  else { return result; }
}

template<class T>
std::vector<T*> Thing::findAll() {
  auto* start = this;  // the starting object will always be the objects from which the method is called
  std::vector<T*> result;
  std::vector<boost::uuids::uuid> scanned; // list of scanned objects by means of UUID
  adder(&scanned,start->getUuid()); // add start to scanned items to avoid looping

  // scan for object of type T in start object relations
  for (auto i : start->relations) {
    auto tested = avoider(start,i);
    adder(&scanned,tested->getUuid()); //add to scanned list
    if (dynamic_cast<T*>(tested)) {
      result.push_back(dynamic_cast<T*>(tested));
    }
  }

  // if not found above, try searching in related objects
  for (auto i : start->relations) {
    auto tested = avoider(start,i);
    adder(&scanned,tested->getUuid()); // add to scanned list
    auto res = looker<T>(tested, &scanned);
    if (res != nullptr) {
      result.push_back(looker<T>(tested, &scanned));
    }
  }

  if (result.size() == 0) {
      std::cout << "Not found." << std::endl;
      abort(); }
  else { return result; }
}

template<class T>
std::vector<T*> Thing::getRelatedObjects() {
   std::vector<T*> objs;

   for(auto i: relations) {
      if(T* r0 = dynamic_cast<T*>(i->getRange())) {
        objs.push_back(r0);
      }
   }

   return objs;
}
