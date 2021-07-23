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

#ifndef PARTICLEBOND_H
#define PARTICLEBOND_H

#include "bond.h"
#include "particle.h"
#include <memory>

/// The most generic bond between particles implementing the methods for the calculation
/// of the sintering level as function of the bond distance and the particle radius.
/// If the particles radius changes for surface condensation, the sintering level is
/// updated accordingly.
/// The bond distance is changed in time only by the Friedlander-Koch sintering law.

template<typename P>
class ParticleBond : public Bond<P> {

    double A; ///< Common area [m2]

public:

    /// Constructor. Initialize the sintering level with a given value.
    /// \param bd initial bond distance
    ParticleBond(std::shared_ptr<P> p0, std::shared_ptr<P> p1, double s=0);

    /// Returns the sintering level calculated by comparing the original bond distance and
    /// the actual particles radii. Sintering level for this class can only change by
    /// particle growth due to surface condensation
    double get_sintering_level() const;

    /// Sintering level update returning the sintering level by means of the
    /// Friedlander-Koch model.
    /// \param dt timestep for the sintering update
    /// \param tau sintering time
    void sintering(double dt, double tau);

    double get_bond_distance() const;

private:

    /// Area of the particle formed by p0 and p1
    /// \param r0 p0 radius [m]
    /// \param r1 p1 radius [m]
    double final_coalescence_area(double r0, double r1) const;

    /// Coefficent for the calculation of the sintering level.
    /// \param r0 p0 radius [m]
    /// \param r1 p1 radius [m]
    double sintering_coeff(double r0, double r1) const;
};

#include "particlebond.cpp"

#endif // PARTICLEBOND_H
