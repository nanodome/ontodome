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

#include "gasphase.h"
#include <iostream>

double GasPhase::get_viscosity() const {
  double visc = 0;
  for (auto i : species) {
      double ci = c[hash.at(i->getName())];
      visc += i->getScalarProperty("viscosity")*ci;
  }
  return visc;
}

double GasPhase::get_p_sat(MolecularEntity* spec, double T) const {

  std::vector<double> p_sat = spec->getVectorProperty("saturation pressure coefficients");

  return 1.01e5 * pow(10.0,(p_sat[0]-(p_sat[1]/T))); ;
}

double GasPhase::get_n_sat(MolecularEntity* spec, double T) const {

  return get_p_sat(spec,T)/(K_BOL*T);
}

double GasPhase::get_S(MolecularEntity* spec, double T) const {

    double m_frac = c[hash.at(spec->getName())];
    double ns = m_frac*p/(K_BOL*T);

    double s_sat = get_n_sat(spec,T);

    return ns/s_sat;
}

double GasPhase::get_s_ten(MolecularEntity* spec, double T) const {

    std::vector<double> s_ten = spec->getVectorProperty("surface tension coefficients");

    return s_ten[0]-s_ten[1]*(T-s_ten[2]);
}

double GasPhase::get_rho() const {

    return get_average_molecular_mass() * p/(K_BOL*T);
}


double GasPhase::get_average_molecular_mass() const {

    double m = 0.0;

    for(size_t i=0; i<species.size(); ++i)
        m += species[i]->getScalarProperty("mass") * c[i];

    return m;
}


double GasPhase::get_gas_flux() const {

    double flux = 0.0;

    for(size_t i=0; i<species.size(); ++i) {

        double m_gas = species[i]->getScalarProperty("mass");

        // n_s * m_s * sqrt(3*K_BOL*T/m_s);
        flux += c[i]*p/(K_BOL*T) * m_gas * sqrt(3*K_BOL*T/m_gas);
    }

    return flux;
}


void GasPhase::print() {

    for(auto& sp : species)
        std::cout << sp->getName() << '\t';
    std::cout << std::endl;

    for(auto& cs : c)
        std::cout << cs << '\t';
    std::cout << std::endl << std::endl;

    std::cout << "T[K]\t" << "p[Pa]\t" << "n[#/m3]" << std::endl;
    std::cout << T << "\t" << p << "\t" << p/(K_BOL*T) << std::endl << std::endl;

}

