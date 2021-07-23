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

#ifndef PARTICLEPHASE_H
#define PARTICLEPHASE_H

#include "nanodome.h"
#include "gasphase/gasphase.h"
#include "nucleationtheory.h"

#include <list>
#include <memory>

/// Concept:
/// A --> Aggregate

template<typename A>
class ParticlePhase {

protected:

    /// Cubic control volume size [m3]
    double volume;

    /// Container for the aggregates
    std::list<std::shared_ptr<A>> aggregates;

public:

    /// Constructor with a volume equivalent to a cube of 1um side
    /// \param _volume volume size [m3]
    ParticlePhase(double _volume) : volume(_volume) {}

    /// Return the size of the control volume that contains the particles [m]
    double get_volume() const { return volume; }

    /// Correct volume due to gas expansion
    /// \param dt timestep [s]
    /// \param gp gas phase surrounding particles
    void volume_expansion(double dt, const GasPhase& gp);

    /// Get the number of aggregates.
    int get_aggregates_number() const { return aggregates.size(); }

    /// Get the aggregates density [#/m3]
    double get_aggregates_density() const { return aggregates.size()/volume; }

    /// Get the mean fractal dimensions of the aggregates.
    double get_mean_fractal_dimension() const;

    /// Get the aggregates mean spherical diameter [m]
    double get_aggregates_mean_spherical_diameter() const;

    /// Get the mean number of particles in an aggregate
    double get_mean_particles_number() const;

    /// smallest particle diameter [m]
    double get_particles_smallest_diameter() const;

	/// biggest particle diameter [m]
	///	\param int _id 
	double get_particles_biggest_diameter(int& _id) const;

    /// Get the mean collision diameter of the aggregates [m]
    double get_mean_collision_diameter() const;

    /// Get the mean sintering level
    double get_mean_sintering_level() const;

    /// Get the particles mean diameter [m]
    double get_particles_mean_diameter() const;

	/// Get all primary particles sizes (for PSD evaluation)
	std::valarray<double> get_particles_sizes() const;

	/// Get all aggregates spherical (for PSD evaluation)
	std::valarray<double> get_aggregates_sizes() const;

protected:

    /// Evaluates condensation for each aggregate and return the volumetric molecules
    /// consumption [#/m3 s]
    /// \param dt timestep [s]
    /// \param Fs surface condensation rate [#/m2 s]
    double condensation(double dt, double Fs);

    /// Execute a sintering step for each aggregate.
    /// \param dt timestep [s]
    /// \param T temperature [K]
    void sintering(double dt, double T);
};

#include "particlephase.cpp"

#endif // PARTICLEPHASE_H
