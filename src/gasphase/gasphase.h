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

#include <vector>
#include <valarray>
#include <map>

#define ND_VISCOSITY 5e-5 // argon viscosity at 1000K


/// Class implementing the gas phase.
/// The gas phase is univocally determined knowing pressure, temperature and species molar fractions.
/// Access to species properties is provided by species name by means of an hash map.
class GasPhase {

protected:

    double p; ///< Gas phase pressure [Pa]
    double T; ///< Gas phase temperature [K]

    double gamma; ///< expansion coefficient [1/s]

    std::map<std::string,std::size_t> hash; ///< Hash map for name-to-index resolution

    // IMPORTANT species and c indexes must be always syncronized
    std::vector<MolecularQuantity*> species; ///< Species container
    std::valarray<double> c; ///< Species molar fraction

public:

    /// Get gas phase pressure [Pa]
    double get_p() const { return p; }

    /// Get single species number density [#/m3]
    double get_n(std::string formula) const { return c[hash.at(formula)]*p/(K_BOL*T); }

	/// Get species molar fraction [%]
	double get_c(std::string formula) const { return c[hash.at(formula)]; }

    /// Get gas phase number density [#/m3]
    double get_n() const { return p/(K_BOL*T); }

    /// Get superaturation ratio
    double get_S(std::string formula) const;

    /// Get the expansion coefficient [1/s]
    double get_gamma() const { return gamma; }

    /// Get gas phase mass density [kg/m3]
    double get_rho() const;

    /// Get gas phase molecules average mass [kg]
    double get_average_molecular_mass() const;

    /// Get species
    Species get_species(std::string formula) const { return species[hash.at(formula)]; }

    /// Get gas phase temperature [K]
    double get_T() const { return T; }

    /// Get gas phase mean free path [m]
    double get_mfp() const { return (ND_VISCOSITY/p) * sqrt(M_PI*K_BOL*T/(2.*get_average_molecular_mass())); }

    /// Get gas flux used for Langevin dynamics [kg/m2 s]
    double get_gas_flux() const;

    /// Print gas phase parameters
    void print();
};

#endif // GASPHASE_H
