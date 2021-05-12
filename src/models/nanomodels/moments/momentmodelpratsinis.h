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

#ifndef MOMENTMODELPRATSINIS_H
#define MOMENTMODELPRATSINIS_H

#include <math.h>
#include <iostream>
#include <fstream>

#include "../../../ontodome.h"
#include "../../../base/thing.h"

/// Implementation of the Pratsinis 1988 model.
/// Simultaneous Nucleation, Condensation and Coagulation in Aerosol Reactors,
/// S.E. Pratsinis, Journal of Colloid and Interface Science, Vol. 124, No. 2, August 1988
class MomentModelPratsinis : public MesoscopicModel {

protected:
    double M0; ///< Nanoparticles number density [#/m3]
    double M1; ///< Nanoparticles total volume [m3/m3]
    double M2; ///< Nanoparticles ... [m6/m3]

    /// attempt to reduce computational time
    double s_mass;
    double s_m_vol;
    double s_bulk_density_liq;
    double s_bulk_density_sol;
    double s_T_melt;

    /// TEST
    double M1_cond;

public:
    /// Standard constructor.
    MomentModelPratsinis() : MesoscopicModel() {
      // List of compatible entities
      createRelationTo<isModelFor,Thing>(new HomonuclearMolecule);
      createRelationTo<isModelFor,Thing>(new HeteronuclearMolecule);

      // List of required models
      createRelationTo<requiresModelFor,Thing>(new GasMixture);
      createRelationTo<requiresModelFor,Thing>(new HomonuclearMolecule);
      createRelationTo<requiresModelFor,Thing>(new HeteronuclearMolecule);

      M0=0.; M1=0.; M2=0.;
    }

    /// Constructor with initialization for the moments
    MomentModelPratsinis(double _M0, double _M1, double _M2) : MesoscopicModel() {
      // List of compatible entities
      createRelationTo<isModelFor,Thing>(new PolyatomicEntity);

      // List of required models
      createRelationTo<requiresModelFor,Thing>(new GasMixture);
      createRelationTo<requiresModelFor,Thing>(new PolyatomicEntity);

      M0 = _M0; M1 = _M1; M2 = _M2;
    }

    /// Initialize the method - attempt to save computational time
    template <class NT, class SPEC> void initialize(NT* nt, SPEC* sp) {

      s_mass = sp -> template getRelatedObject<Mass>()[0] -> template getRelatedObject<Scalar>()[0]->data;
      s_m_vol = nt -> template get_m_volume();
      s_bulk_density_liq = sp -> template getRelatedObject<BulkDensityLiquid>()[0] -> template getRelatedObject<Scalar>()[0] -> data;
      s_bulk_density_sol = sp -> template getRelatedObject<BulkDensitySolid>()[0] -> template getRelatedObject<Scalar>()[0] -> data;
      s_T_melt = sp -> template getRelatedObject<MeltingPoint>()[0] -> template getRelatedObject<Scalar>()[0] -> data;
    }

    /// Nanoparticle density [#/m3]
    double get_n_density() const { return M0; }

    /// Nanoparticle mean diameter [m]
    double get_mean_diameter() const { return (M0>0.0) ? pow(6*M1/(M_PI*M0),1./3.) : 0.0; }

    /// Nanoparticle volume density [m3/m3]
    double get_total_volume() const { return M1; }

    /// ... [m6/m3]
    double get_total_area() const { return M2; }

    /// Timestep calculation
    /// \param dt timestep size [s]
    /// \param T temperature [K]
    /// \param J nucleation rate [#/m3 s]
    /// \param j stable cluster size [#]
    /// \param S supersaturation ratio
    template <class GM, class NT, class SPEC> double timestep(double dt,  GM* gm, NT* nt, SPEC* sp) {

      double T = gm-> template get_T();
      double S = gm -> template get_S(sp, T);
      double gamma = gm -> template get_gamma();
      double t_csi1 = csi1(sp,gm,T);
      double t_csi2 = csi2(T);

      /*if (S < 1.0)
              S = 1.0;*/

      double J = nt -> template nucleation_rate(sp, gm, T);
      double j = nt -> template stable_cluster_size(sp, gm, T);

      // temporary moments values
      double M0_, M1_, M2_  = 0;

      // volume of the critical size cluster
      // if S<1 the NucleationTheory object should return 1 as stable critical cluster size
      double vm = j * s_m_vol;

      /// dissipation limit
      double vm1 = 0.0;

      // nucleation source for each moment
      double M0_nucl = J;
      double M1_nucl = J * vm;
      double M2_nucl = J * vm*vm;

      double vg = geometric_mean_v();
      double ln2sg = ln2_standard_dev();

      // condensation source for each moment
      double M1_cond = 0;
      double M2_cond = 0;;
      M1_cond = t_csi1*(S - 1)*M_k(2. / 3., ln2sg, vg);
      M2_cond = 2 * t_csi1*(S - 1)*M_k(5. / 3., ln2sg, vg);

      // TEST
//      M1_cond = 0.0;
//      M2_cond = 0.0;

      double sg = exp(sqrt(ln2sg));

      // coagulation source for each moment
      double M0_coag =   t_csi2*zeta0(sg)*(  M_k(2./3.,ln2sg,vg) * M_k(-1./2.,ln2sg,vg) +
                                            2*M_k(1./3.,ln2sg,vg) * M_k(-1./6.,ln2sg,vg) +
                                              M_k(1./6.,ln2sg,vg) * M0);

      double M2_coag = 2*t_csi2*zeta2(sg)*(  M_k(5./3.,ln2sg,vg) * M_k(1./2.,ln2sg,vg) +
                                            2*M_k(4./3.,ln2sg,vg) * M_k(5./6.,ln2sg,vg) +
                                              M_k(7./6.,ln2sg,vg) * M1);

      // TEST
      M0_coag = 0.0;
      M2_coag = 0.0;

      // dissolution source: if S<1 then the flux of particles becoming smaller than a two monomers cluster
      // due to evaporation is removed from the particle set
      // in this case vm = 1 must hold
      double sigma = 3.0 * sqrt(ln2sg);
      double M0_diss = 0;
      double M1_diss = 0;
      double M2_diss = 0;

      if((S<1)&&(vg>0)&&(sigma>0)) {

        //vm *= 2;
        //
        //double nm = M0/(sigma*sqrt(2*M_PI)) * exp(-pow(log(vm/vg),2)/(2*sigma*sigma)) / vm;
        //      double G0 = csi1(T)*pow(vm,2./3.)*(S-1);

        //      M0_diss = -G0*nm; // particle number reduction
        //      M1_diss = -G0*nm*vm; // particle volume reduction
        //      M2_diss = -G0*nm*vm*vm; // M2 reduction

        vm1 = vg*std::exp(-1.82*1.41*sigma);

        double nm = M0 / (sigma*sqrt(2 * M_PI)) * exp(-pow(log(vm1 / vg), 2) / (2 * sigma*sigma)) / vm1;
        double G0 = t_csi1*pow(vm1, 2. / 3.)*(S - 1);

        M0_diss = -G0*nm; // particle number reduction
        M1_diss = -G0*nm*vm1; // particle volume reduction
        M2_diss = -G0*nm*vm1*vm1; // M2 reduction

        //M0_diss *= 0;
        //M1_diss *= 1.0e5;
        //M2_diss *= 1.0e5;

      }

      // moments equations solution (simple first order explicit method) (TO BE IMPROVED!!!)
      M0_ = M0 + (M0_nucl           - M0_coag - M0_diss) * dt - gamma*M0*dt;
      M1_ = M1 + (M1_nucl + M1_cond           - M1_diss) * dt - gamma*M1*dt;
      M2_ = M2 + (M2_nucl + M2_cond + M2_coag - M2_diss) * dt - gamma*M2*dt;

//      M0_ = M0 + (M0_nucl - M0_coag - M0_diss) * dt;
//      M1_ = M1 + (M1_nucl + M1_cond - M1_diss) * dt;
//      M2_ = M2 + (M2_nucl + M2_cond + M2_coag - M2_diss) * dt;

      // nucleating species consumption for
      // + homogeneous nucleation: J*j
      // + heterogenous nucleation: M1_cond/vm
      // - dissolution of evaporating particles reducing their size below vm
      double g = J*j + M1_cond/s_m_vol - M1_diss/s_m_vol;

      // TEST
      //g = 0.0;

      // DEBUG
      if (std::isnan(M0_) || std::isinf(M0_)) {
              std::cout
                      << " M0: " << M0 << std::endl
                      << " M1: " << M1 << std::endl
                      << " M2: " << M2 << std::endl
                      << " M0_nucl: " << M0_nucl << std::endl
                      << " M0_coag: " << M0_coag << std::endl
                      << " M0_diss: " << M0_diss << std::endl
                      << " csi2(T): " << t_csi2 << std::endl
                      << " zeta0(sg): " << zeta0(sg) << std::endl
                      << " sg: " << sg << std::endl
                      << " vg: " << vg << std::endl
                      << " ln2sg: " << ln2sg << std::endl
                      << " M_k(2. / 3., ln2sg, vg): " << M_k(2. / 3., ln2sg, vg) << std::endl
                      << " M_k(-1. / 2., ln2sg, vg): " << M_k(-1. / 2., ln2sg, vg) << std::endl
                      << " M_k(1. / 3., ln2sg, vg): " << M_k(1. / 3., ln2sg, vg) << std::endl
                      << " M_k(-1. / 6., ln2sg, vg): " << M_k(-1. / 6., ln2sg, vg) << std::endl
                      << " M_k(1. / 6., ln2sg, vg): " << M_k(1. / 6., ln2sg, vg) << std::endl
                      << std::endl;
              std::cout
          << "J: " << J << std::endl;
//              system("PAUSE"); //does not work on Linux-based systems
              std::cin.get();
        }

       // TEST
      if (M0_<0.0 || M1_ < 0.0 || M2_ < 0.0){
              M1_ = 0.0;
              M2_ = 0.0;
              M0_ = 0.0;
      }

      // update moments
      M0 = M0_; M1 = M1_; M2 = M2_;

      return g;
    }

    /// Return the geometric mean volume [m3]
    double geometric_mean_v() {

      return ((M0>0.0)&&(M2>0.0)) ? M1*M1/(pow(M0,1.5)*sqrt(M2)) : 0.0;
    }

    /// Return the standard deviation of the lognormal distribution
    double ln2_standard_dev() {

        double value = 0.0;

        double arg = M0*M2/(M1*M1);

        if(arg>1.0) {
            value = (1./9.) * log(arg);
            // Debug
            if (std::isinf(value)) {
              std::cout << "SG COMPUTATION" << std::endl;
              std::cout
                  << " M0: " << M0 << std::endl
                  << " M1: " << M1 << std::endl
                  << " M2: " << M2 << std::endl;
            }
        }

        return value;
    }

    /// Return the geometric standard deviation of the lognormal distribution
    double standard_dev() { return 3.0 * sqrt(ln2_standard_dev()); }

    /// Return the condenstion source term
    template <class SPEC, class NT> double get_cond_term(SPEC* species, NT* nt) const { return M1_cond / nt-> template get_m_volume(species); }

    /// Get M0
    double get_M0()const { return M0; }

    /// Get M1
    double get_M1()const { return M1; }

    /// Get M2
    double get_M2()const { return M2; }

    /// Get M_k
    double M_k(double k, double ln2sg, double vg) {

        return (vg>0.0) ? M0*pow(vg,k)*exp(4.5*k*k*ln2sg) : 0.0;
    }

    double zeta0(double sg) {

        return 0.633 + 0.092*sg*sg - 0.022*sg*sg*sg;
    }

    double zeta2(double sg) {

        return 0.39 + 0.5*sg - 0.214*sg*sg + 0.029*sg*sg*sg;
    }


    template <class SPEC, class GM> double csi1(SPEC* species, GM* gm, double T) {
      double s_n_sat = gm -> template get_n_sat(species,T);

      return s_m_vol*s_n_sat*4.835975862049408*sqrt(K_BOL*T/(2*M_PI*s_mass));
    }

    double csi2(double T) {
      double s_bulk = get_bulk_density(T);

      return 0.787623317899743*sqrt(6*K_BOL*T/s_bulk);
    }

    /// Get Species bulk density [kg/m3]
    /// \param T Species temperature [K]
    double get_bulk_density(double T) {
      return (T<s_T_melt) ? s_bulk_density_sol : s_bulk_density_liq;
    }

    /// Get Lognormal
    void get_lognormal_val(const std::string& path) {

      std::ofstream _o_file;
      _o_file.open(path, std::ofstream::app);

      // Moment 0
      (std::isnan(M0) || std::isinf(M0)) ? _o_file << 0.0 : _o_file << M0;
      _o_file << " ";

      // Moment 1
      (std::isnan(M1) || std::isinf(M1)) ? _o_file << 0.0 : _o_file << M1;
      _o_file << " ";

      // Moment 2
      (std::isnan(M2) || std::isinf(M2)) ? _o_file << 0.0 : _o_file << M2;
      _o_file << " ";

      //vg
      double vg = geometric_mean_v();
      (std::isnan(vg) || std::isinf(vg)) ? _o_file << 0.0 : _o_file << vg;
      _o_file << " ";

      /// sigma g
      double ln2sg = ln2_standard_dev();
      (std::isnan(ln2sg) || std::isinf(ln2sg)) ? _o_file << 0.0 : _o_file << ln2sg;
      _o_file << " ";


      _o_file << std::endl;

      // close file
      _o_file.close();

    }

    /// Get Lognormal (Streamlines)
    void get_lognormal_val(int _s_index) {

      std::ofstream _o_file;
      _o_file.open("Lognormal_"+std::to_string(_s_index) + ".dat", std::ofstream::app);

      // Moment 0
      (std::isnan(M0) || std::isinf(M0)) ? _o_file << 0.0 : _o_file << M0;
      _o_file << " ";

      // Moment 1
      (std::isnan(M1) || std::isinf(M1)) ? _o_file << 0.0 : _o_file << M1;
      _o_file << " ";

      // Moment 2
      (std::isnan(M2) || std::isinf(M2)) ? _o_file << 0.0 : _o_file << M2;
      _o_file << " ";

      //vg
      double vg = geometric_mean_v();
      (std::isnan(vg) || std::isinf(vg)) ? _o_file << 0.0 : _o_file << vg;
      _o_file << " ";

      /// sigma g
      double ln2sg = ln2_standard_dev();
      (std::isnan(ln2sg) || std::isinf(ln2sg)) ? _o_file << 0.0 : _o_file << ln2sg;
      _o_file << " ";


      _o_file << std::endl;

      // close file
      _o_file.close();
    }

};


#endif // MOMENTMODELPRATSINIS_H