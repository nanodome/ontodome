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

#ifndef AGGREGATE_H
#define AGGREGATE_H

#include "mesoobject.h"
#include "objectcounter.h"
#include "particlebond.h"

#include <list>
#include <memory>
#include <valarray>


/// The Aggregate class is a container of connected particles.
/// It provides methods for particles accretion due to surface condensation.
/// Sintering level calculation is provided by the bonds.
/// This type provides also a method for coalescence if sintering level of a bond is
/// beyond a specific limit (e.g. s > 0.95).
///
/// Concepts:
/// P --> Particle

template<typename P>
class Aggregate : public MesoObject, public ObjectCounter<Aggregate<P>> {

protected:

    /// Particles composing the aggregate
    std::list<std::shared_ptr<P>> particles;

    /// Bonds connecting the particles
    std::list<std::shared_ptr<ParticleBond<P>>> bonds;

public:

    /// Constructor
    /// \param p0 pointer to the newly created particle.
    Aggregate(std::shared_ptr<P> p0) { particles.push_back(p0); }

	/// Constructor from a list of particles and bonds
	Aggregate(std::list<std::shared_ptr<P>> _particles,
		std::list<std::shared_ptr<ParticleBond<P>>> _bonds) {
		particles = _particles;
		bonds = _bonds;
	}

    /// Copy constructor
    Aggregate(const Aggregate& a1);

    /// Get aggregate mass as sum of the particles masses [kg]
    double get_mass() const final;

    /// Get aggregate volume as sum of the particles volumes [m3]
    double get_volume() const final;

    /// Number of particles in the aggregate [#]
    double get_particles_number() const { return particles.size(); }

    /// Get the smallest particle diameter in the aggregate [m]
    double get_particles_smallest_diameter() const;

	/// Get the biggest particle diameter in the aggregate [m]
	double get_particles_biggest_diameter() const;

    /// Get the number of bonds in this aggregate
    int get_bonds_number() const { return bonds.size(); }

    /// Particles mean diameter [m]
    double get_particles_mean_diameter() const;

    /// Aggregate area as sum of the particles surface area [m2]
    double get_surface_area() const;

    /// Get the equivalent spherical surface of the aggregate [m2]
    double get_spherical_surface() const;

    /// Get the diameter of the equivalent spherical surface [m]
    double get_spherical_diameter() const;

    /// Mean sintering level of the aggregate
    double get_mean_sintering_level() const;

    /// Get particles number accounting for the reduction in size due to sintering
    double get_reduced_particles_number() const;

    /// Aggregate active area taking into account reduction due to sintering [m2]
    double get_active_surface_area() const;

	/// Get primary particles diameter
	std::valarray<double> get_particles_diameter();

    /// Surface condensation on the particles of this aggregate. (Condensation is
    /// calculated using full particle surface, and is not scaled using the active surface
    /// of the aggregate - TODO)
    /// \param dt timestep size [s]
    /// \param Fs surface condensation rate [#/m2 s]
    double condensation(double dt, double Fs);

    /// Execute a sintering step for each particle in the aggregate.
    /// \param dt timestep [s]
    /// \param T temperature [K]
    void sintering(double dt, double T);

    /// Check if a bond has a sintering level higher than a threshold and provide
    /// coalescence for the particles.
    void coalescence(double max_sint_level);

	///	Merges two aggregates
	///	\param std::shared_ptr<Aggregate<P>> a1 aggregate to merge
	///	\param std::shared_ptr<P> p0 contact particle 1
	///	\param std::shared_ptr<P> p1 contact particle 2
    void merge(std::shared_ptr<Aggregate<P>> _a1,
               std::shared_ptr<P> _p0, std::shared_ptr<P> _p1);

	///	Merges two aggregates in a periodic domain
	///	\param std::shared_ptr<Aggregate<P>> a1 aggregate to merge
	///	\param std::valarray<double> shift coordinates shift
	///	\param std::shared_ptr<P> p0 contact particle 1
	///	\param std::shared_ptr<P> p1 contact particle 2
	void merge_periodic(std::shared_ptr<Aggregate<P>> a1,
						std::valarray<double> shift,
						std::shared_ptr<P> _p0,
						std::shared_ptr<P> _p1);

    void merge(std::shared_ptr<Aggregate<P>> a1);
};

#include "aggregate.cpp"

#endif // AGGREGATE_H
