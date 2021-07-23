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

#include "pbmaggregate.h"


template<typename P>
double PBMAggregate<P>::get_collision_diameter() const {

    double d_avg = this->get_particles_mean_diameter();
    double n = this->get_particles_number();

    return d_avg*pow(n,1.0/D_f);
}

template<typename P>
double PBMAggregate<P>::get_n_monomers() const{

	double tot_n = 0.0;
    for (auto p = this->particles.begin(); p != this->particles.end(); p++) {
		tot_n += (*p)->get_n();
	}
	return tot_n;

}

template<typename P>
Species PBMAggregate<P>::get_species() const{

	
     return (*this->particles.begin())->get_species();

}
