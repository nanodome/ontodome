#ifndef SPECIES_H
#define SPECIES_H

#include <string>
#include <cmath>
#include <valarray>
#include <boost/algorithm/cxx11/any_of.hpp>

#include "../base/object.h"

class Species : public Object
{
  // molecular properties
/*  String formula(std::string form); ///< Substance formula (e.g. "ZnO")
  String name(std::string name);  */  ///< IUPAC name (e.g. "zinc monoxide")

//  ScalarQuantity mass(double mass, std::string unit); ///< molecular mass [kg]
//  ScalarQuantity molecular_volume(double molvol, std::string unit); ///< molecular volume [m3]

//  ScalarQuantity sigma(double sigma, std::string unit); ///< L-J sigma value [m]
//  ScalarQuantity eps(double eps, std::string unit); ///< L-J epsilon value [J]

//  // substance properties
//  ScalarQuantity density(double dens, std::string unit); ///< Density [kg/m3]
//  ScalarQuantity T_melt(double T_melt, std::string unit); ///< Melting temperature [K]
//  ScalarQuantity viscosity(double visc, std::string unit); ///< viscosity at 300K [kg/m s]

//  // surface tension is expressed as s_ten_A - s_ten_B*(T - s_ten_C) [N/m]
//  ScalarQuantity s_ten_A(double s_ten_A, std::string unit); ///< Surface tension coefficient A [N/m]
//  ScalarQuantity s_ten_B(double s_ten_B, std::string unit); ///< Surface tension coefficient B [N/(m K)]
//  ScalarQuantity s_ten_C(double s_ten_C, std::string unit); ///< Surface tension coefficient C [K]

//  // saturation pressure is expressed as log10(psat) = (psat_A - psat_B/T) * 1.01e5 [Pa]
//  ScalarQuantity p_sat_A(double p_sat_A, std::string unit); ///< Saturation pressure coefficient A
//  ScalarQuantity p_sat_B(double p_sat_B, std::string unit); ///< Saturation pressure coefficient B [K]

//  // NASA polynomials for thermodynamic properties
//  ScalarQuantity nasa_min_t(double nasa_min_t, std::string unit); ///< Min temperature range [K]
//  ScalarQuantity nasa_max_t(double nasa_max_t, std::string unit); ///< Max temperature range [K]
//  ScalarQuantity nasa_mid(double nasa_mid, std::string unit);   ///< Common temperature between hi and low range [K]
//  VectorQuantity nasa_ht(std::valarray<double> nasa_ht, std::string unit); ///< High temperature range polynomials
//  VectorQuantity nasa_lt(std::valarray<double> nasa_lt, std::string unit); ///< Low temperature range polynomials

  public:
    Species();

//    template<typename T> auto getProperty(this, std::string name){
//      if (!boost::algorithm::any_of_equal(this->properties->first, name)){
//        return _prop;
//        }
//      return _prop;
//    }

    void addProperty(double _prop, std::string _unit, std::string _name);
};

#endif // SPECIES_H
