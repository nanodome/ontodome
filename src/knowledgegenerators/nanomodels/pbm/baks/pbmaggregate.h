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

#ifndef PBMAGGREGATE_H
#define PBMAGGREGATE_H

#include "fractalaggregate.h"

/// An aggregate type for which the fractal dimension is prescribed upon creation.
///
/// Concept:
/// P --> Particle

template<typename P>
class PBMAggregate : public FractalAggregate<P> {

    /// Fractal dimension of the aggregate.
    double D_f;

public:

    /// Constructor. The fractal dimension is needed as initial assumption for
    /// the aggregate.
    /// \param _D_f fractal dimension of the aggregate
    /// \param p0 pointer to the first particle composing the aggregate
    PBMAggregate(double _D_f, std::shared_ptr<P> p0)
        : FractalAggregate<P>(p0), D_f(_D_f) {}

    /// Aggregate fractal dimension
    double get_fractal_dimension() const { return D_f; }

    /// Get the collision diameter [m]
    double get_collision_diameter() const;

	/// Get number of monomers composing the aggregate [MONOSPECIES]
	double get_n_monomers() const;

	/// Get aggregates
    std::list<std::shared_ptr<P>> get_particles() const { return this->particles; }

	/// Get the species composing the aggregate [MONOSOECIES]
	Species get_species()const;
};

#include "pbmaggregate.cpp"

#endif // PBMAGGREGATE_H
