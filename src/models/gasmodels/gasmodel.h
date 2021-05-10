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

#ifndef GASMODEL_H
#define GASMODEL_H

#include <valarray>

#include "../../ontodome.h"
#include "../../base/thing.h"

/// Class implementing the gas phase.
/// The gas phase is univocally determined knowing pressure, temperature and species molar fractions.
/// Access to species properties is provided by hasPart relations defined in the GasMoxture object
class GasModel  : public ContinuumModel {

protected:

    double p; ///< Gas phase pressure [Pa]
    double dpdt; ///< Gas phase pressure time derivative [Pa/s]
    double T; ///< Gas phase temperature [K]
    double dTdt; ///< Gas phase temperature time derivative [K/s]
    double gamma; ///< expansion coefficient [1/s]
    std::valarray<double> w; ///< Species molar fractions

public:
    GasModel()
    {
      createRelationTo<isModelFor,Thing>(new GasMixture);
      createRelationTo<isModelFor,Thing>(new HomonuclearMolecule);
      createRelationTo<isModelFor,Thing>(new HeteronuclearMolecule);
    }

    std::string getClassName() const { return "Gas Continuum Model"; }

    /// Run the model
    void run(GasMixture* gm, double dt, std::valarray<double> w_cons) {
      // Get the GasMixture conditions
      p = get<Pressure,GasMixture>(gm);
      dpdt = get<PressureTimeDerivative,GasMixture>(gm);
      T = get<Temperature,GasMixture>(gm);
      dTdt = get<TemperatureTimeDerivative,GasMixture>(gm);
      w = get_molar_fractions(gm);

      // Check if molar fraction vector total sum is smaller than 1 for consistency
      if (w.sum() > 1.) { abort(); }

      // Check if consumption vector and available species number matches
      // If not, the code must abort
      std::vector<PolyatomicEntity*> specs = gm->getRelatedObject<PolyatomicEntity>();
      if (specs.size() != w_cons.size()) { abort(); }

      // Compute the new species number densities
      double n = p/(K_BOL*T); // number of monomers on m3
      double w_cons_tot = w_cons.sum(); // total number of consumed monomers

      std::valarray<double> ns = w*n; // number of monomers for each species

      gamma = w_cons_tot/n + dTdt/T - dpdt/p;

      // simple explicit ODE timestep
      // equation is solved for the number density
      for (std::size_t i = 0; i < w_cons.size(); ++i) {
          ns[i] += (w_cons[i] - ns[i] * gamma) * dt; // Volume contraction
      }

      // gas phase pressure and temperature update
      T += dTdt*dt;
      p += dpdt*dt;

      // update the GasMixture
      update<Pressure, GasMixture>(gm,p);
      update<Temperature, GasMixture>(gm,T);

      // molar fractions update
      n = p/(K_BOL*T);
      w = ns/n;

      // update molar fractions of each species
      for (std::size_t i = 0; i < w.size(); ++i) {
          update<MolarFraction, PolyatomicEntity>(specs[i],w[i]);
      }

      print(gm);
    }

    /// Get gas phase molar fractions
    std::valarray<double> get_molar_fractions(GasMixture* gm)
    {
      std::vector<PolyatomicEntity*> specs = gm->getRelatedObject<PolyatomicEntity>();

      if (!specs.empty())
      {
          w.resize(specs.size());
          for (std::size_t i = 0; i < specs.size(); ++i)
          {
            w[i] = specs[i]->getRelatedObject<MolarFraction>()[0]->getRelatedObject<Scalar>()[0]->data;
          }
          return w;
      }
      else { abort(); }
    }

    /// Get gas phase properties
    template <class T, class OBJ> double get(OBJ* obj) const
    {
      double val;
      auto vals = obj-> template getRelatedObject<T>();
      if (!vals.empty())
      {
          if (vals.size() > 1) {
              abort();
          }
          else
          {
              val = vals[0]->template getRelatedObject<Scalar>()[0]->data;
          }
          return val;
      }
      else { abort(); }
    }

    /// Update gas phase properties
    template <class T, class OBJ> void update(OBJ* obj, double val) const
    {
      auto vals = obj->template getRelatedObject<T>();
      if (!vals.empty())
      {
          if (vals.size() > 1) {
              abort();
          }
          else
          {
              vals[0]->template getRelatedObject<Scalar>()[0]->data = val;
          }
      }
      else { abort(); }
    }

    /// Get the expansion coefficient [1/s]
    double get_gamma() const { return gamma; }

    /// Get gas phase molecules average viscosity [Pa s]
    double get_average_viscosity(GasMixture* gm) const {
      double visc = 0;

      std::vector<PolyatomicEntity*> specs = gm->getRelatedObject<PolyatomicEntity>();
      for(size_t i=0; i<w.size(); ++i)
       visc += specs[i]->getRelatedObject<Viscosity>()[0]->getRelatedObject<Scalar>()[0]->data * w[i];

      return visc;
    }

    /// Get gas phase number density [#/m3]
    double get_n() const { return p/(K_BOL*T); }

    /// Get superaturation ratio [#]
    /// \param spec selected species
    /// \param T temperature [K]
    template <class TT> double get_S(TT* spec, double T) const {

      double m_frac = spec->template getRelatedObject<MolarFraction>()[0]->template getRelatedObject<Scalar>()[0]->data;
      double ns = m_frac*p/(K_BOL*T);

      double s_sat = get_n_sat<TT>(spec,T);
      return ns/s_sat;
    }

    /// Saturation density [#/m3]
    /// \param spec selected species
    /// \param T temperature [K]
    template <class TT> double get_n_sat(TT* spec, double T) const { return spec->template getRelatedObject<SaturationPressure>()[0]->get_p_sat(T) /(K_BOL*T); }

    /// Get gas phase mass density [kg/m3]
    double get_density(GasMixture* gm) const { return get_average_molecular_mass(gm) * p/(K_BOL*T); }

    /// Get gas phase molecules average mass [kg]
    double get_average_molecular_mass(GasMixture* gm) const {

      double m = 0.;

      std::vector<PolyatomicEntity*> specs = gm->getRelatedObject<PolyatomicEntity>();
      for(size_t i=0; i<w.size(); ++i)
          m += specs[i]->getRelatedObject<Mass>()[0]->getRelatedObject<Scalar>()[0]->data * w[i];

      return m;
    }

    /// Get gas phase mean free path [m]
    double get_mfp(GasMixture* gm) const { return (get_average_viscosity(gm)/p) * sqrt(M_PI*K_BOL*T/(2.*get_average_molecular_mass(gm))); }

    /// Get gas flux used for Langevin dynamics [kg/m2 s]
    double get_gas_flux(GasMixture* gm) const {
      double flux = 0.;

      std::vector<PolyatomicEntity*> specs = gm->getRelatedObject<PolyatomicEntity>();
      for(size_t i=0; i<w.size(); ++i) {

          double m_gas = specs[i]->getRelatedObject<Mass>()[0]->getRelatedObject<Scalar>()[0]->data;

          // n_s * m_s * sqrt(3*K_BOL*T/m_s);
          flux += w[i]*p/(K_BOL*T) * m_gas * sqrt(3*K_BOL*T/m_gas);
      }

      return flux;
    }

    /// Print gas phase parameters
    void print(GasMixture* gm) {
      for(auto& sp : gm->getRelatedObject<PolyatomicEntity>())
          std::cout << sp->getRelatedObject<IUPAC>()[0]->data << '\t';
      std::cout << "Sum" << '\t';
      std::cout << std::endl;

      for(auto& cs : w)
          std::cout << cs << '\t';
      std::cout << w.sum() << '\t';
      std::cout << std::endl << std::endl;

      std::cout << "T[K]\t" << "p[Pa]\t" << "n[#/m3]" << std::endl;
      std::cout << T << "\t" << p << "\t" << p/(K_BOL*T) << std::endl << std::endl;
    };

};

#endif // GASMODEL_H
