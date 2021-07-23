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

#include "particlebond.h"
#include <cmath>


template<typename P>
ParticleBond<P>::ParticleBond(std::shared_ptr<P> p0, std::shared_ptr<P> p1, double s)
    : Bond<P>(p0,p1) {

    double r0 = 0.5*this->get_v0()->get_diameter();
    double r1 = 0.5*this->get_v1()->get_diameter();

    double Af = final_coalescence_area(r0,r1);
    double coeff = sintering_coeff(r0,r1);

    // calculate area depending on sintering level
    A = Af/(coeff + (1.0 - coeff)*s);
}


template<typename P>
double ParticleBond<P>::get_sintering_level() const {

    double r0 = 0.5*this->get_v0()->get_diameter();
    double r1 = 0.5*this->get_v1()->get_diameter();

    double Af = this->final_coalescence_area(r0,r1);

    // calculate sintering level
    double coeff = this->sintering_coeff(r0,r1);

    // return sintering level
    return (Af/A - coeff)/(1.0 - coeff);
}


template<typename P>
void ParticleBond<P>::sintering(double dt, double tau) {

    double r0 = 0.5*this->get_v0()->get_diameter();
    double r1 = 0.5*this->get_v1()->get_diameter();

    double Af = final_coalescence_area(r0,r1);

    // simple explicit timestep
    A += (Af-A)*(dt/tau);

    // check for consistency in the result
    if(A<Af) A = Af;
}

template<typename P>
double ParticleBond<P>::get_bond_distance() const {

    // linear interpolation between minimum and maximum sintering distance

    double r0 = 0.5*this->get_v0()->get_diameter();
    double r1 = 0.5*this->get_v1()->get_diameter();

    double dmax = r0+r1;
    double dmin = fabs(r1-r0);

    double s = get_sintering_level();
	if (s < 0.0)
		s = 0.0; 
	
	return (1- s)*dmax + s*dmin;
}


template<typename P>
double ParticleBond<P>::final_coalescence_area(double r0, double r1) const {

    // volume of the final particle
    double Vf = 4./3. * M_PI * (r0*r0*r0 + r1*r1*r1);

    // surface of the final coalesced particle
    return 1.464591887561523*pow(6.*Vf,2./3.); // pow(M_PI,1./3.)
}


template<typename P>
double ParticleBond<P>::sintering_coeff(double r0, double r1) const {

    double alpha = r0/r1;

    // this coefficient is suited for particles with different diameter,
    // for particles with same diameter would be coeff = pow(2.0,-1./3.)
    return pow(1+alpha*alpha*alpha,2./3.)/(1+alpha*alpha);
}
