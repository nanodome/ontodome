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

#ifndef CNT_H
#define CNT_H

#include <cmath>
#include <limits>
#include <exception>
#include <stdexcept>
#include <iostream>

#include "../../../ontodome.h"
#include "../../../base/thing.h"

/// Classical Nucleation Theory (CNT) implementation. This class is the implementation of the
/// CNT whose rates are calculated using the properties of the condensing species.
class ClassicalNucleationTheory : public MesoscopicModel {

public:

    /// Default constructor.
    /// \param _species species type
    ClassicalNucleationTheory() {
      // List of compatible entities
      createRelationTo<isModelFor,Thing>(new HomonuclearMolecule);
      createRelationTo<isModelFor,Thing>(new HeteronuclearMolecule);

      // List of required models
      createRelationTo<requiresModelFor,Thing>(new GasMixture);
    };

    /// Primary particles formation rate [#/m3 s]
    /// \param T temperature [K]
    /// \param S supersaturation ratio
    template <class TT, class GasM> double nucleation_rate(TT* species, GasM* gasmodel, double T) const {

      double rate = 0.;
      //unsigned int fp_control_state = _controlfp(_EM_INEXACT, _MCW_EM);

      // Get the supersaturation ration from GasPhase
      double S = gasmodel-> template get_S(species,T);

      // check if species is saturated, if not nucleation rate is left to zero
      if(S>1.0) {

          // normalized surface tension
          double theta = species-> template getRelatedObject<SurfaceTension>()[0]-> template get_s_ten(T)*get_m_surface(species)/(K_BOL*T);
          double ns_sat = gasmodel-> template get_n_sat(species,T);

          double A = (S*ns_sat*ns_sat);
          double B = get_m_volume(species);
          double C1 = 2.0*species-> template getRelatedObject<SurfaceTension>()[0]->get_s_ten(T) ;
          double C2 = M_PI*species-> template getRelatedObject<Mass>()[0]->template getRelatedObject<Scalar>()[0]->data;
          double D1 = theta - 4.0*pow(theta, 3) / (27.0*pow(log(S), 2));
          double C = sqrt(C1 / C2);
          double D = exp(D1);

          rate = A * B * C * D;

      } else if (S==1) {

          rate = std::numeric_limits<double>::min();
      }

      return rate;
    }

    /// Smallest stable cluster particle number [#]
    /// \param T temperature [K]
    /// \param S supersaturation ratio
    template <class TT, class GasM> double stable_cluster_size(TT* species, GasM* gasmodel, double T) const {

      double c_size = 1.0;

      // Get the supersaturation ration from GasPhase
      double S = gasmodel-> template get_S(species,T);

      // check if species is saturated; if not, the stable size has no sense and is set to one
      // meaning that the smallest cluster is a single monomer (no-cluster)
      if(S>1) {
          c_size = 2.0 * get_m_surface(species) * species-> template getRelatedObject<SurfaceTension>()[0]->get_s_ten(T) / (3*K_BOL*T*log(S));
          c_size = pow(c_size,3);
                  }

      return c_size;
    }

    /// Smallest stable cluster diameter [m]
    /// \param ns nucleating species concentration [#/m3]
    /// \param T temperature [K]
    /// \param S supersaturation ratio
    template <class TT, class GasM> double stable_cluster_diameter(TT* species, GasM* gasmodel, double T) const {

      // Get the supersaturation ration from GasPhase
      double S = gasmodel-> template get_S(species,T);

      return 2*pow( (3./4.) * get_m_volume(species) * stable_cluster_size(species,gasmodel,S) / M_PI, 1./3.);
    }

    /// Surface condensation rate [#/m2 s]
    /// \param T temperature [K]
    /// \param S supersaturation ratio
    template <class TT, class GasM> double condensation_rate(TT* species, GasM* gasmodel, double T) const{

      // Get the supersaturation ration from GasPhase
      double S = gasmodel-> template get_S(species,T);

      return species-> template getRelatedObject<SaturationPressure>()[0]->get_p_sat(T)*(S-1.0) / sqrt(2*M_PI*species-> template getRelatedObject<Mass>()[0]->template getRelatedObject<Scalar>()[0]->data*K_BOL*T);
    }

    /// Molecular Volume [m3]
    /// \param Selected species
    template <class TT> double get_m_volume(TT* species) const {
      return species-> template getRelatedObject<Mass>()[0]->template getRelatedObject<Scalar>()[0]->data/(AMU*N_AVO*1000*species-> template getRelatedObject<BulkDensityLiquid>()[0]->template getRelatedObject<Scalar>()[0]->data);
    }

    /// Molecular Surface [m2]
    /// \param Selected species
    template <class TT> double get_m_surface(TT* species) const { return 4.8360*pow(get_m_volume(species),2./3.); }

};

#endif // CNT_H
