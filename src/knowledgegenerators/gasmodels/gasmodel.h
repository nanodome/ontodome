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
#include <map>

#include "../../ontodome.h"

/// Class implementing the gas phase.
/// The gas phase is univocally determined knowing pressure, temperature and species molar fractions.
/// Access to species properties is provided by hasPart relations defined in the GasMixture object.
class GasModel  : public SoftwareModel {

protected:
    double p; ///< Gas phase pressure [Pa].
    double T; ///< Gas phase temperature [K].
    double dpdt; ///< Gas phase pressure time derivative [Pa/s].
    double dTdt; ///< Gas phase temperature time derivative [K/s].
    double gamma; ///< expansion coefficient [1/s].
    std::vector<SingleComponentComposition*> specs; ///< vector containing all the species.
    std::valarray<double> w; ///< Species molar fractions.
    std::map<std::string,std::size_t> hash; ///< Hash map for name-to-index resolution.
    bool init = false; ///< Boolean value which stores whether the model is initialized or not.

public:

    GasModel()
    {
//      createRelationTo<hasModel,GasMixture>(new GasMixture);
    }

    std::string getClassName() const { return "Gas Continuum Model"; }

    /// Initialize the model by looking for the required inputs through the relations graph.
    void initialize() {
      // Get the GasMixture conditions
      auto gp = this->findNearest<GasMixture>();
      specs = gp->getRelatedObjects<SingleComponentComposition>();
      p = gp->findNearest<Pressure>()->findNearest<Real>()->data;
      dpdt = gp->findNearest<PressureTimeDerivative>()->findNearest<Real>()->data;
      T = gp->findNearest<Temperature>()->findNearest<Real>()->data;
      dTdt = gp->findNearest<TemperatureTimeDerivative>()->findNearest<Real>()->data;
      w = get_molar_fractions();

      for(size_t i=0; i<specs.size(); ++i)
          hash[specs[i]->name] = i;

      // Mark the model as initialized
      init = true;
    }

    /// Run the model for a given timestep and species consumption array.
    /// \param dt temporal step [s].
    /// \param w_cons species consumption array [#/m3/s].
    void timestep(double dt, std::valarray<double> w_cons) {
      // Initialize the model if not done before
      if (init == false) { initialize(); }

      // Check if molar fraction vector total sum is smaller than 1 for consistency
      // A 1% excess is tolerated
      if (w.sum() >= 1.*1.01) { abort(); }

      // Check if consumption vector and available species number matches
      // If not, the code must abort
      if (specs.size() != w_cons.size()) { abort(); }

      // Compute the new species number densities
      double n = get_n(); // number of monomers on m3
      double w_cons_tot = w_cons.sum(); // total number of consumed monomers

      std::valarray<double> ns = w*n; // number of monomers for each species

      gamma = w_cons_tot/n + dTdt/T - dpdt/p;

      // simple explicit ODE timestep
      // equation is solved for the number density
      for (std::size_t i = 0; i < w_cons.size(); ++i) {
          ns[i] += (w_cons[i] - ns[i] * gamma) * dt; // Volume contraction
          if (ns[i] < 0) abort();
      }

      // gas phase pressure and temperature update
      T += dTdt*dt;
      p += dpdt*dt;

      // molar fractions update
      n = get_n();
      w = ns/n;
    }

//    /// Push final state to GasMixture
//    void add_temporal_state(Matter* gp,double time) {
//      //create the state to push
//      auto state = new Matter;

//      //add the current GasMixture thermodynamic state
//      state->createRelationsTo<hasProperty,Quantity>({
//        new Pressure (new Scalar(p), new Unit("Pa")),
//        new Temperature (new Scalar(T), new Unit("K")),
//        new PressureTimeDerivative (new Scalar(dpdt), new Unit("Pa/s")),
//        new TemperatureTimeDerivative (new Scalar(dTdt), new Unit("K/s"))});

//      //update and add the current species state
//      for (auto i : specs) {
//          update<MolarFraction,Matter>(i,w[hash.at(i->getUuid())]);
//      }
//      state->createRelationsTo<hasPart,PolyatomicEntity>(specs);

//      //push the state to GasMixture
//      gp->push_state(state,time);
//    }

    /// Get the gas phase molar fractions array.
    std::valarray<double> get_molar_fractions()
    {
      if (!specs.empty())
      {
          w.resize(specs.size());
          for (std::size_t i = 0; i < specs.size(); ++i)
          {
            w[i] = specs[i]->mol->findNearest<Real>()->data;
          }
          return w;
      }
      else { abort(); }
    }

//    /// Update gas phase properties
//    template <class T, class OBJ> void update(OBJ* obj, double val) const
//    {
//      auto vals = obj->template getRelatedObjects<T>();
//      if (!vals.empty())
//      {
//          if (vals.size() > 1) {
//              abort();
//          }
//          else
//          {
//              vals[0]->template getRelatedObjects<Real>()[0]->data = val;
//          }
//      }
//      else { abort(); }
//    }

    /// Get the expansion coefficient [1/s].
    double get_gamma() const { return gamma; }

    /// Get the temperature [K].
    double get_T() const { return T; }

    /// Get the pressure [pa].
    double get_p() const { return p; }

    /// Get gas phase molecules average viscosity [Pa s].
    double get_average_viscosity() const {
      double visc = 0;

      for(size_t i=0; i<w.size(); ++i)
       visc += specs[i]->getRelatedObjects<Viscosity>()[0]->getRelatedObjects<Real>()[0]->data * w[i];

      return visc;
    }

    /// Get gas phase number density [#/m3].
    double get_n() const { return p/(K_BOL*T); }

    /// Get the superaturation ratio [#].
    /// \param spec selected species.
    /// \param T temperature [K].
    template <class TT> double get_S(TT* spec, double T) const {

      double m_frac = w[hash.at(spec->name)];
      double ns = m_frac*get_n();
      double n_sat = get_n_sat<TT>(spec,T);

      return ns/n_sat;
    }

    /// Get the saturation density [#/m3].
    /// \param spec selected species.
    /// \param T temperature [K].
    template <class TT> double get_n_sat(TT* spec, double T) const {
      spec->template getRelatedObjects<SaturationPressure>()[0]->template getRelatedObjects<SaturationPressureMaterialRelation>()[0]->run();
      return spec->template getRelatedObjects<SaturationPressure>()[0]->template getRelatedObjects<Real>()[0]->data /(K_BOL*T); }

    /// Get the gas phase mass density [kg/m3].
    double get_density() const { return get_average_molecular_mass() * p/(K_BOL*T); }

    /// Get the gas phase molecules average mass [kg].
    double get_average_molecular_mass() const {

      double m = 0.;

      for(size_t i=0; i<w.size(); ++i)
          m += specs[i]->getRelatedObjects<Mass>()[0]->getRelatedObjects<Real>()[0]->data * w[i];

      return m;
    }

    /// Get the gas phase mean free path [m].
    double get_mfp() const { return (get_average_viscosity()/p) * sqrt(M_PI*K_BOL*T/(2.*get_average_molecular_mass())); }

    /// Get the gas flux used for Langevin dynamics [kg/m2 s].
    double get_gas_flux() const {
      double flux = 0.;

      for(size_t i=0; i<w.size(); ++i) {

          double m_gas = specs[i]->getRelatedObjects<Mass>()[0]->getRelatedObjects<Real>()[0]->data;

          // n_s * m_s * sqrt(3*K_BOL*T/m_s);
          flux += w[i]*p/(K_BOL*T) * m_gas * sqrt(3*K_BOL*T/m_gas);
      }

      return flux;
    }

    /// Prints the gas phase most relevant properties.
    void print() {
      for(auto& sp : specs)
          std::cout << sp->name << '\t';
      std::cout << "Sum" << '\t';
      std::cout << std::endl;

      for(auto& cs : w)
          std::cout << cs << '\t';
      std::cout << w.sum() << '\t';
      std::cout << std::endl << std::endl;

      std::cout << "T[K]\t" << "p[Pa]\t" << "n[#/m3]\t  " << "gamma[1/s]\t" << std::endl;
      std::cout << T << "\t" << p << "\t" << get_n() << "\t" << "  " << gamma << "\t" << std::endl << std::endl;
    };

};

#endif // GASMODEL_H
