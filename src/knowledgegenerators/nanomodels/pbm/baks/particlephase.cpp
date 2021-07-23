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

#include "particlephase.h"
#include "ndm_random.h"

#include<cfloat>


template<typename A>
void ParticlePhase<A>::volume_expansion(double dt, const GasPhase& gp) {

        volume *= exp(gp.get_gamma()*dt);
}


template<typename A>
double ParticlePhase<A>::get_mean_fractal_dimension() const {

    if(aggregates.size()==0)
        return 0;

    double Df = 0;

    for(auto& a: aggregates)
        Df += a->get_fractal_dimension();

    return Df/aggregates.size();
}


template<typename A>
double ParticlePhase<A>::get_aggregates_mean_spherical_diameter() const {

    if(aggregates.size()==0)
        return 0;

    double d = 0;

    for(auto& a: aggregates)
        d += a->get_spherical_diameter();

    return d/aggregates.size();
}

template<typename A>
double ParticlePhase<A>::get_mean_particles_number() const {

    if(aggregates.size()==0)
        return 0;

    double n = 0;

    for(auto& a: aggregates)
        n += a->get_particles_number();

    return n/aggregates.size();
}


template<typename A>
double ParticlePhase<A>::get_particles_smallest_diameter() const {

    double d;
    double d_min = 1e20;

    for(auto& a: aggregates) {
        d = a->get_particles_smallest_diameter();
        if(d<d_min) d_min = d;
    }

    return d_min;
}

template<typename A>
double ParticlePhase<A>::get_particles_biggest_diameter(int& _id) const {

	double d;
	double d_max = DBL_MIN;
	_id = -1;

	for (auto& a : aggregates) {
		d = a->get_particles_biggest_diameter();
		if (d > d_max) {
			d_max = d; _id = a->get_id();
		}
	}

	return d_max;
}

template<typename A>
double ParticlePhase<A>::get_mean_collision_diameter() const {

    if(aggregates.size()==0)
        return 0;

    double d = 0;

    for(auto& a: aggregates)
        d += a->get_collision_diameter();

    return d/aggregates.size();
}


template<typename A>
double ParticlePhase<A>::get_mean_sintering_level() const {

    if(aggregates.size()==0)
        return 0;

    double s = 0;

    for(auto& a: aggregates)
        s += a->get_mean_sintering_level();

    return s/aggregates.size();
}


template<typename A>
double ParticlePhase<A>::condensation(double dt, double Fs) {

    double n_tot = 0;

	/*for (auto& a : aggregates) {
		n_tot += a->condensation(dt, Fs);
	}*/

	for (auto a = this->aggregates.begin(); a != this->aggregates.end(); ++a) {
		n_tot += (*a)->condensation(dt, Fs);

		// if, due negative condensation, the aggregate disappear, erase it from the list
		if ((*a)->get_particles_number() == 0)
			a = this->aggregates.erase(a);

		if (a == this->aggregates.end())
			break;
	}
	
	return n_tot/(volume*dt);
}


template<typename A>
void ParticlePhase<A>::sintering(double dt, double T) {

    for(auto& a: aggregates)
        a->sintering(dt,T);
}

template<typename A>
double ParticlePhase<A>::get_particles_mean_diameter() const {

    double n = 0;
    for (auto& a : this->aggregates)
        n += a->get_particles_number();

    double d = 0;
    for (auto& a : this->aggregates) {
        for (auto& p : a->get_particles_diameter()) {
            if (p == 0){
                n -= 1;
                continue;}
            else{
                d += p;}
        }
    }

    if (n == 0)
        return 0;
    else
        return d/n;
}

template<typename A>
std::valarray<double> ParticlePhase<A>::get_particles_sizes() const {

    double n = 0;

	for (auto& a : this->aggregates)
		n += a->get_particles_number();

    std::valarray<double> sizes(n);

    int idx = 0;
	for (auto& a : this->aggregates) {
		for (auto& p : a->get_particles_diameter()) {
            sizes[idx] = p;
			idx++;
		}
	}

    return sizes;
}

template<typename A>
std::valarray<double> ParticlePhase<A>::get_aggregates_sizes() const {

	int n = aggregates.size();
	std::valarray<double> sizes(n);

	int idx = 0;
	for (auto& a : this->aggregates) {
		sizes[idx] = a->get_spherical_diameter();
		idx++;
	}

	return sizes;
}
