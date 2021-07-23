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

#ifndef PBMPARTICLEPHASE_H
#define PBMPARTICLEPHASE_H

#include "particlephase.h"
#include "../aggregate/pbmaggregate.h"

/// Concept:
/// A --> PBMAggregate

template<typename A>
class PBMParticlePhase : public ParticlePhase<A> {

public:

    /// Maximum timestep size [s]
    double dt_max;

    /// Maximum number of aggregates. Further nucleations will lead to the removal of
    /// an aggregate and a volume reduction.
    int max_aggregates_number;

    /// Minumum number of aggregates. Further reductions will lead to
    /// the duplication of a randomly chosen aggregate.
    int min_aggregates_number;

    /// The particle phase is initialized by defining its volume.
    /// \param _D_f Fractal dimension
    /// \param _volume PP volume [m3]
    PBMParticlePhase(double _volume = 1e-18) : ParticlePhase<A>(_volume) {}

    /// Get total particle processes rate [#/s]
    double get_total_processes_rate(const GasPhase& gp, const NucleationTheory& nt, const Species& sp);

    /// Perform the PP timestep [#/m3 s]
    double timestep(double dt, const GasPhase& gp, const NucleationTheory& nt, const Species& sp);

	/// Get the aggregate with index n
	///	\param int index [0, N particles - 1]
	std::shared_ptr<PBMAggregate<Particle>> get_aggregate(int index);

	/// Get all particles diameters
	std::vector<double> get_all_particles_diameters();

	/// Get all aggregates diameters
	std::vector<double> get_all_aggregates_diameters();

	/// Sets maximum number of aggregates
	void set_max_aggregates(int _max) { this->max_aggregates_number = _max; }

	/// Sets minimum number of aggregates
	void set_min_aggregates(int _min) { this->min_aggregates_number = _min; }

	/// Calculate the coagulation rate using the pre-calculated parameters [#/s]
	double coagulation_rate(double T) const;

	//  parameters and functions for majorant kernel calculations
	double K1, K2, K21, K22;

	// DUMP all aggregates to a file
	void dump_aggregates(std::ofstream& _dump);

private:

    // particles process rates [1/s]
    double R_nucl, R_coag, R_null, R_tot;

    /// Adjust the aggregates number between N_max and N_min.
    void aggregates_number_balance();

    /// Recalculate the rates for nucleation and coagulation.
    void update_rates(const GasPhase& gp, const NucleationTheory& nt, const Species& sp);

    /// Nucleation of a single new aggregate [#/m3 s]
    virtual double nucleation(double j, Species s) = 0;

    /// Coagulation of two randomly chosen aggregates
    void coagulation();

// ##########################################################
////  parameters and functions for majorant kernel calculations
//    double K1,K2,K21,K22;
//
    /// Updates K1, K2, K21, K22 using actual aggregates population.
    void update_majorant_kernel_parameters();

    

    /// Select two aggregates for coagulation.
    void aggregates_selection(std::shared_ptr<A>& a0, std::shared_ptr<A>& a1);
// ###########################################################
};

#include "pbmparticlephase.cpp"

#endif // PBMPARTICLEPHASE_H
