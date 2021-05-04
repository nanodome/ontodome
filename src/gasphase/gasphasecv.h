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

#ifndef GASPHASECV_H
#define GASPHASECV_H

#include "gasphase.h"

class GasPhaseCV : public GasPhase {

public:

    /// Standard constructor.
    /// \param _p pressure [Pa]
    /// \param _T temperature [K]
    /// \param _species string vector with species names
    /// \param _c vector with molar fraction of each species ordered as _species
    GasPhaseCV(double _p, double _T, std::vector<Species> _species, std::valarray<double> _c);

    /// Explicit constant pressure solver
    /// \param dt timestep [s]
    /// \param dTdt temperature gradient [K/s]
    /// \param dpdt pressure gradient [Pa/s]
    /// \param w vector with species source terms [#/m3 s]
    void timestep(double dt, double dTdt, double dpdt, std::valarray<double> w = {0.0,0.0});
};

#endif // GASPHASECV_H
