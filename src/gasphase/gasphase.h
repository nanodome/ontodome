/*
    NanoDome - H2020 European Project NanoDome, GA n.646121
    (www.nanodome.eu, https://github.com/nanodome/nanodome)
    e-mail: Emanuele Ghedini, emanuele.ghedini@unibo.it

    Copyright (C) 2018  Alma Mater Studiorum - Università di Bologna

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GASPHASE_H
#define GASPHASE_H

#include "../base/thing.h"
#include "../nanodome.h"

#include <vector>
#include <valarray>
#include <map>

/// Class implementing the gas phase.
/// The gas phase is univocally determined knowing pressure, temperature and species molar fractions.
/// Access to species properties is provided by species name by means of an hash map.
class GasPhase  : public ContinuumModel
{

protected:

    double p; ///< Gas phase pressure [Pa]
    double T; ///< Gas phase temperature [K]
    double gamma; ///< expansion coefficient [1/s]

    std::map<std::string,std::size_t> hash; ///< Hash map for name-to-index resolution

    // IMPORTANT species and c indexes must be always syncronized
    std::vector<MolecularEntity*> species; ///< Species container
    std::valarray<double> c; ///< Species molar fraction

public:
    GasPhase() : ContinuumModel("GasPhase") {}

    std::string getClassName() const { return "Gas Phase Continuum Model"; }

    /// Get gas phase pressure [Pa]
    double get_p() const { return p; }

    /// Get gas phase temperature [K]
    double get_T() const { return T; }

    /// Get the expansion coefficient [1/s]
    double get_gamma() const { return gamma; }

    /// Get the expansion coefficient [1/s]
    double get_viscosity() const;

    /// Get single species number density [#/m3]
    double get_n(std::string name) const { return c[hash.at(name)]*p/(K_BOL*T); }

    /// Get species molar fraction [%]
    double get_c(std::string name) const { return c[hash.at(name)]; }

    /// Get gas phase number density [#/m3]
    double get_n() const { return p/(K_BOL*T); }

    /// Get superaturation ratio
    double get_S(MolecularEntity* spec, double T) const;

    /// Saturation pressure [Pa]
    /// \param T temperature [K]
    double get_p_sat(MolecularEntity* spec, double T) const;

    /// Saturation density [#/m3]
    /// \param T temperature [K]
    double get_n_sat(MolecularEntity* spec, double T) const;

    /// Get species surface tension [N/m]
    double get_s_ten(MolecularEntity* spec, double T) const;

    /// Get gas phase mass density [kg/m3]
    double get_rho() const;

    /// Get gas phase molecules average mass [kg]
    double get_average_molecular_mass() const;

    /// Get species
    MolecularEntity* get_species(std::string name) const {

      MolecularEntity* spec;
      for (auto i : species) {
          if ( i->getName() == name) {
              spec = i;
          } else {
              std::cout << "The species " << name << "is missing. Please add it." << std::endl;
              abort();
          }
      }
      return spec;
    }

    /// Get gas phase mean free path [m]
    double get_mfp() const { return (get_viscosity()/p) * sqrt(M_PI*K_BOL*T/(2.*get_average_molecular_mass())); }

    /// Get gas flux used for Langevin dynamics [kg/m2 s]
    double get_gas_flux() const;

    /// Print gas phase parameters
    void print();
};

#endif // GASPHASE_H
