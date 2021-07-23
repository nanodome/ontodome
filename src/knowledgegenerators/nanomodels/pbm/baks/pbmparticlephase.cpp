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

#include "pbmparticlephase.h"
#include "particle.h"
#include "nanodome.h"
#include "ndm_random.h"

#include <iostream>


template<typename A>
double PBMParticlePhase<A>::get_total_processes_rate(const GasPhase& gp, const NucleationTheory& nt, const Species& sp) {

    // setting the simulation parameters with default values
    dt_max = 1e-8; // [s]
	//max_aggregates_number = 2000;
    //min_aggregates_number = 1990;

    update_rates(gp,nt,sp);

    return R_tot;
}

template<typename A>
std::shared_ptr<PBMAggregate<Particle>> PBMParticlePhase<A>::get_aggregate(int index) {

	std::list<std::shared_ptr<PBMAggregate<Particle>>>::iterator it;
	if (index >= this->aggregates.size()) {
		//error
		//return NULL;
	}
	else{
		it = this->aggregates.begin();
		std::advance(it, index);
	}

	return (*it);
}


template<typename A>
double PBMParticlePhase<A>::timestep(double dt, const GasPhase& gp, const NucleationTheory& nt, const Species& sp) {

    // MONOSPECIES!!!

    // nucleating species concentration and temperature
    double T  = gp.get_T();
    double S  = gp.get_S(sp.get_formula());

    // monomers consuption rate [#/m3 s]
    double consumption = 0;

    // choose between nucleation and coagulation processes
    double rho = ndm::uniform_double_distr(ndm::rand_gen);

    if(rho <= R_nucl/R_tot) {

        // nucleation
        double j = nt.stable_cluster_size(T,S);

        consumption += nucleation(j,sp)/dt;

    } else if(rho<=(R_nucl+R_coag)/R_tot) {

        coagulation();
    }

    // apply condensation
	if (S < 1) S = 1.0;
    double Fs = nt.condensation_rate(T,S);
	
    consumption += this->condensation(dt,Fs);

    // apply sintering
    this->sintering(dt,T);

    // adjust the aggregates number
    aggregates_number_balance();

    return consumption;
}


template<typename A>
void PBMParticlePhase<A>::coagulation() {

    // choose two aggregates a0!=a1
    std::shared_ptr<A> a0,a1;
    do { aggregates_selection(a0,a1); } while (a0==a1);

    // calculate kernels
    double m0 = a0->get_mass();
    double m1 = a1->get_mass();
    double d0 = a0->get_collision_diameter();
    double d1 = a1->get_collision_diameter();

    // free molecular kernel
    double K = sqrt(1/m0 + 1/m1) * square(d0 + d1);

    // free molecular majorant kernel
    double k_maj = 2.0;
    double K_maj = k_maj * (1/sqrt(m0) + 1/sqrt(m1)) * (square(d0) + square(d1));

    // check if it is a real or a fictious coagulation
    double delta = ndm::uniform_double_distr(ndm::rand_gen);
    if(delta <= K/K_maj)
        a0->merge(a1);

    this->aggregates.remove(a1);
}


template<typename A>
void PBMParticlePhase<A>::update_rates(const GasPhase& gp, const NucleationTheory& nt, const Species& sp) {

    // nucleating species concentration and temperature
    double T  = gp.get_T();
    double S  = gp.get_S(sp.get_formula());

    // nucleation parameters for this timestep
    double J = nt.nucleation_rate(T,S);

    R_coag = R_null = 0;
	double volume_tmp = this->get_volume();
	R_nucl = J * volume_tmp;

    // calculate the coagulation reaction rates
    if(this->aggregates.size() >= 2) {
        update_majorant_kernel_parameters();
        R_coag = coagulation_rate(T);
    }

    // total reaction rates
    R_tot = R_coag + R_nucl;

    // if particle process rate is lower than a specified 1/DT_MAX then provide an
    // artificial NULL process to increase the timestep length
    if(R_tot < 1.0/dt_max) {
        R_null = 1.0/dt_max - R_tot;
        R_tot += R_null;
    }
}


template<typename A>
void PBMParticlePhase<A>::update_majorant_kernel_parameters() {

    double d_col;
    double mass;

    K1 = K2 = K21 = K22 = 0;

    for(auto& a: this->aggregates) {

        d_col = a->get_collision_diameter();
        mass = a->get_mass();

        d_col = square(d_col);
        mass = 1./sqrt(mass);

        K1 += d_col*mass;
        K21 += d_col;
        K22 += mass;
    }

    K1 *= 2.0 * (this->aggregates.size() - 2.0);
    K2 = 2.0 * K21 * K22;


}


template<typename A>
double PBMParticlePhase<A>::coagulation_rate(double T) const {

    double k_maj = 2;

    return k_maj * sqrt(0.5*M_PI*K_BOL*T)*(K1+K2)/this->get_volume();
}


template<typename A>
void PBMParticlePhase<A>::aggregates_selection(std::shared_ptr<A>& a0, std::shared_ptr<A>& a1) {

    double rho = ndm::uniform_double_distr(ndm::rand_gen);

    if(rho<K1/(K1+K2)) {

        // pick a random a0 with w = 1
        auto it = this->aggregates.begin();

        double w0_rnd = ndm::uniform_double_distr(ndm::rand_gen);

        a0 = *(std::next(it,int(w0_rnd*(this->aggregates.size()-1))));

        // pick a random a1 with w = d^2/sqrt(m) using Cumulative Density Funct
        double w1 = 0;
        double w1_rnd = ndm::uniform_double_distr(ndm::rand_gen);
        double w1_norm = 0.5*K1/(this->aggregates.size()-2.0);

        for(auto it = this->aggregates.begin(); it!=this->aggregates.end(); ++it) {

            double d_col = (*it)->get_collision_diameter();
            double mass = (*it)->get_mass();

            d_col = square(d_col); mass = 1./sqrt(mass);
            w1 += d_col*mass/w1_norm;

            if(w1>w1_rnd) {
                a1 = *(it);
                break;
            }
        }

    } else {

        double w0 = 0;
        double w1 = 0;

        double w0_rnd = ndm::uniform_double_distr(ndm::rand_gen);
        double w1_rnd = ndm::uniform_double_distr(ndm::rand_gen);

        double w0_norm = K21;
        double w1_norm = K22;

        // pick a random a0 with w = d^2
        for(auto it=this->aggregates.begin(); it!=this->aggregates.end(); ++it) {

            double d_col = (*it)->get_collision_diameter();

            d_col = square(d_col);
            w0 += d_col/w0_norm;

            if(w0>w0_rnd) {
                a0 = *(it);
                break;
            }
        }

        // pick a random a1 with w = 1/sqrt(m)
        for(auto it=this->aggregates.begin(); it!=this->aggregates.end(); ++it) {

            double mass = (*it)->get_mass();

            mass = 1./sqrt(mass);
            w1 += mass/w1_norm;

            if(w1>w1_rnd) {
                a1 = *(it);
                break;
            }
        }
    }
}

// DUMP all aggregates to a file
template<typename A>
void PBMParticlePhase<A>::dump_aggregates(std::ofstream& _dump) {

	_dump << "agg_id" << '\t' << "agg_mass" << '\t' 
		<< "agg_vol" << '\t' << "col_dim" << '\t' 
		<< "mean_diam" << '\t' << "particles num"
		<< "DF"
		<< std::endl;
	for (auto agg = this->aggregates.begin(); agg != this->aggregates.end(); agg++) {
		_dump << (*agg)->get_id() 
			<< " " << (*agg)->get_mass()
			<< " " << (*agg)->get_volume()
			<< " " << (*agg)->get_collision_diameter()
			<< " " << (*agg)->get_particles_mean_diameter()
			<< " " << (*agg)->get_particles_number()
			<< " " << (*agg)->get_fractal_dimension()
			<< " " << std::endl;
	}

	for (auto agg = this->aggregates.begin(); agg != this->aggregates.end(); agg++) {
		auto particles = (*agg)->get_particles();
		_dump << "ID: " << (*agg)->get_id() << std::endl;
		_dump << "particles: " << std::endl;
		for (auto p = particles.begin(); p != particles.end(); p++) {
			_dump << " mass: " << (*p)->get_mass() << " vol: " << (*p)->get_volume() << std::endl;
		}
		_dump << std::endl;
	}

}

template<typename A>
void PBMParticlePhase<A>::aggregates_number_balance() {

    int N = this->aggregates.size();

    if(N>max_aggregates_number) {

        auto it = this->aggregates.begin();
        double rho = ndm::uniform_double_distr(ndm::rand_gen);

        std::advance(it,int(rho*(N-1)));

        this->aggregates.erase(it);

        this->volume *= (N-1.0)/N;
    }

    if(N==min_aggregates_number) {

        while(N!=max_aggregates_number) {

            auto it = this->aggregates.begin();
            double rho = ndm::uniform_double_distr(ndm::rand_gen);

            std::advance(it,int(rho*(N-1)));
            this->aggregates.emplace_back(new A(**it));

            this->volume *= (N+1.0)/N;
            N++;
        }
    }
}


template<typename A>
std::vector<double> PBMParticlePhase<A>::get_all_particles_diameters() {

	std::vector<double> res;
	
    for (auto agg = this->aggregates.begin(); agg != this->aggregates.end(); agg++) {
		auto particles = (*agg)->get_particles();
		for (auto p = particles.begin(); p != particles.end(); p++) {
			res.push_back((*p)->get_diameter());
		}
	}

	return res;
}

template<typename A>
std::vector<double> PBMParticlePhase<A>::get_all_aggregates_diameters() {

	std::vector<double> res;

    for (auto agg = this->aggregates.begin(); agg != this->aggregates.end(); agg++) {
		res.push_back((*agg)->get_spherical_diameter());
	}

	return res;

}
