#ifndef GASPHASE_H
#define GASPHASE_H

#include "../base/object.h"

#include <vector>
#include <valarray>
#include <map>
#include <iostream>

/// Class implementing the gas phase.
/// The gas phase is univocally determined knowing pressure, temperature and species molar fractions.
/// Access to species properties is provided by species name by means of an hash map.
class GasPhase
{
  protected:
    double p;
    double T;

    double gamma;

    std::map<std::string,std::size_t> hash;

    std::valarray<Object*> species;
    std::valarray<double> n;

  public:
    GasPhase();

//    for (auto i: species){
//          std::cout << i << std::endl;
//    }
};

#endif // GASPHASE_H
