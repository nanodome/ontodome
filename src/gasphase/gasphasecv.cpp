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

#include "gasphasecv.h"


GasPhaseCV::GasPhaseCV(double _p, double _T, std::vector<Species> _species, std::valarray<double> _c) {

    gamma = 1.0;

    p = _p;
    T = _T;
    species = _species;
    c = _c;

    for(size_t i=0; i<species.size(); ++i)
        hash[species[i].get_formula()] = i;
}


void GasPhaseCV::timestep(double dt, double dTdt, double dpdt, std::valarray<double> w) {

    double n = p/(K_BOL*T); // number of monomers on m3
    double wtot = w.sum(); // total number of consumed monomers

    std::valarray<double> ns = c*n; // number of monomers for each species

    gamma = wtot/n + dTdt/T - dpdt/p;

    // simple explicit ODE timestep
    // equation is solved for the number density
	for (std::size_t i = 0; i < c.size(); ++i) {
		ns[i] += (w[i] - ns[i] * gamma) * dt; // Volume contraction
		//ns[i] += (w[i]) * dt; // constant volume
	}

    // gas phase pressure and temperature update
    T += dTdt*dt;
    p += dpdt*dt;

    // molar concentration update
    n = p/(K_BOL*T);
    c = ns/n;
}
