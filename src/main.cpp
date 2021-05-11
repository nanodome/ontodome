#include <iostream>
#include <vector>

#include "ontodome.h"
#include "base/thing.h"
#include "models/gasmodels/gasmodel.h"
#include "models/gasmodels/gasmodelcv.h"
#include "models/nanomodels/nucleation/cnt.h"
#include "models/nanomodels/moments/momentmodelpratsinis.h"

void MomentsRun(GasModel* gm, GasMixture* gp, ClassicalNucleationTheory* cnt, MomentModelPratsinis* mm) {

  // Condensing species datas
  auto spec = gp->getRelatedObject<HomonuclearMolecule>()[0];
  auto spec_name = spec->getRelatedObject<IUPAC>()[0]->data;

  // Initialize the models
  gm->initialize(gp);
  cnt->initialize(spec);
  mm->initialize(cnt,spec);

  // simulation settings
  const double dt = 1e-8; //5e-7;
  double t = 0.;
  double t_end = 0.5;
  int iter = 0;

  const int PRINT_EVERY = 1000;

//  const int SAVE_EVERY = 1000;

//  // Clean lognormal Ag
//  std::string lognormal_path = savep + "Lognormal" +  ".dat";
//  std::ofstream lognormal_f;
//  lognormal_f.open(lognormal_path);
//  lognormal_f.clear();
//  lognormal_f.close();

//  // Plot file Ag
//  std::ofstream plot_data;
//  std::string filename = savep + "MOMENTS_plot" + ".dat";
//  plot_data.open(filename);
//  // Print headlines
//  plot_data << "Time[sec]" << '\t'
//            << "Temp[K]" << '\t'
//            << "SuperSaturation ratio" << '\t'
//            << "Nucl Rate" << '\t'
//            << "Species # density" << '\t'
//            << "Stable cluster size[m]" << '\t'
//            << "AVG Part Num[#]: " << '\t'
//            << "Sint level[%]" << '\t'
//            << "AVG diameter[m]" << '\t'
//            << "Agg. #[#]" << '\t'
//            << "Agg density[#/m3]" << '\t'
//            << "Volume[m]" << '\t'
//            << "AVG fract dim" << '\t'
//            << "M1 cond. value" << '\t'
//            << "M1" << '\t'
//            << "M2" << std::endl;
//    //        << "ts exec time"
//    //        << std::endl;

  // loop over timesteps until the final temperature
  while(t <= t_end) {

      double T = gp->getRelatedObject<Temperature>()[0]->getRelatedObject<Scalar>()[0]->data;
      double S = gm->get_S(spec,T);
//      double J = cnt->nucleation_rate(spec,gm,T);
//      double j = cnt->stable_cluster_size(spec,gm,T);

//      if (T < 300.) { break; }
      if (T <= 300.) {
          gm->update<TemperatureTimeDerivative, GasMixture>(gp,0.);
      }

      // moment method timestep and species consumption retrieval
      double g_cons = mm->timestep(dt, gm, cnt, spec);

      // updating the gas phase
      gm->timestep(gp,dt,{-g_cons,0.});

      ++iter;
      t+=dt;

      if (counter_trigger(iter, PRINT_EVERY)) {
          gm->print();

          std::cout
              << "Time= "<< t << '\t'
              << "Temp= " << T << '\t'
              << "Sat_" << spec_name << "= " << S << '\t'
              << "w_" << spec_name << "= " << gm->get_molar_fractions()[0] << '\t'
              << spec_name << "_cons= " << -g_cons << '\t'
              << "Mean Diam_" << spec_name << "= " << mm->get_mean_diameter() << '\t'
              << "M0_" << spec_name << "= " << mm->get_M0() << '\t'
              << "M1_" << spec_name << "= " << mm->get_M1() << '\t'
              << "M2_" << spec_name << "= " << mm->get_M2() << '\t'
              << std::endl << std::endl;

//            std::cout << "Iteration number: " << iter << std::endl; //per capire quante volte stampa
      }

//      // File savings
//      if (counter_trigger(iter, SAVE_EVERY)) {

//          // save lognormal values
//          mm.get_lognormal_val(lognormal_path);
//          // save simulation data
//          plot_data
//              << t << '\t'
//              << T << '\t'
//              << S << '\t'
//              << J << '\t'
//              << ns << '\t'
//              << j << '\t'
//              << 0.0 << '\t' //dummy particles number
//              << 0.0 << '\t' // dummy sintering level
//              << mm.get_mean_diameter() << '\t'
//              << 0.0 << '\t' // dummy aggregates number
//              << mm.get_density() << '\t'
//              << 0.0 << '\t' // dummy Volume
//              << 0.0 << '\t' // dummy fractal dimension
//              << 0.0 <<  '\t' // dummy interval time
//              << mm.get_cond_term() << '\t'
//              << mm.get_M1() << '\t'
//              << mm.get_M2()
//              << std::endl;
//      }
  }

//  plot_data.close();
}

int main()
{
    GasMixture gp;
    Pressure p(new Scalar(101325.), new Unit("Pa"));
    Temperature T(new Scalar(4000.), new Unit("K"));
    PressureTimeDerivative dpdt(new Scalar(0.), new Unit("Pa/s"));
    TemperatureTimeDerivative dTdt(new Scalar(-1e+5), new Unit("K/s"));

    double ww = 0.0005;
    HomonuclearMolecule Si;
    MolarFraction msi(new Scalar(ww), new Unit("#"));
    IUPAC si("Si");
    Mass masi(new Scalar(28.085*AMU), new Unit("kg"));
    Viscosity musi(new Scalar(7e-5), new Unit("Pa s"));
    BulkDensityLiquid bdlsi(new Scalar(2570.), new Unit("kg/m3"));
    BulkDensitySolid bdssi(new Scalar(2329.), new Unit("kg/m3"));
    MeltingPoint mpsi(new Scalar(1687.), new Unit("K"));
    SaturationPressure psatsi(new Vector({7.5341,23399.}),new Unit("Pa"));
    SurfaceTension stensi(new Vector({0.732,0.000086,1685.}),new Unit("N/m"));
    Si.createRelationTo<hasProperty,MolarFraction>(&msi);
    Si.createRelationTo<hasProperty,IUPAC>(&si);
    Si.createRelationTo<hasProperty,Mass>(&masi);
    Si.createRelationTo<hasProperty,Viscosity>(&musi);
    Si.createRelationTo<hasProperty,BulkDensityLiquid>(&bdlsi);
    Si.createRelationTo<hasProperty,BulkDensitySolid>(&bdssi);
    Si.createRelationTo<hasProperty,MeltingPoint>(&mpsi);
    Si.createRelationTo<hasProperty,SaturationPressure>(&psatsi);
    Si.createRelationTo<hasProperty,SurfaceTension>(&stensi);

//    HeteronuclearMolecule CH4;
//    MolarFraction mch4(new Scalar(0.05), new Unit("#"));
//    IUPAC ch4("CH4");
//    Mass mach4(new Scalar(180*AMU), new Unit("kg"));
//    Viscosity much4(new Scalar(5e-5), new Unit("Pa s"));
//    SaturationPressure psatch4(new Vector({0,0}),new Unit("#"));
//    SurfaceTension stench4(new Vector({0.,0.,0.}),new Unit("#"));
//    CH4.createRelationTo<hasProperty,MolarFraction>(&mch4);
//    CH4.createRelationTo<hasProperPart,IUPAC>(&ch4);
//    CH4.createRelationTo<hasProperty,Mass>(&mach4);
//    CH4.createRelationTo<hasProperty,Viscosity>(&much4);
//    CH4.createRelationTo<hasProperty,SaturationPressure>(&psatch4);
//    CH4.createRelationTo<hasProperty,SurfaceTension>(&stench4);

    HomonuclearMolecule He;
    MolarFraction mhe(new Scalar(1-ww), new Unit("#"));
    IUPAC he("He");
    Mass mahe(new Scalar(4.002602*AMU), new Unit("kg"));
    Viscosity muhe(new Scalar(5e-5), new Unit("Pa s"));
    SaturationPressure psathe(new Vector({0,0}),new Unit("#"));
    SurfaceTension stenhe(new Vector({0.,0.,0.}),new Unit("#"));
    He.createRelationTo<hasProperty,MolarFraction>(&mhe);
    He.createRelationTo<hasProperty,IUPAC>(&he);
    He.createRelationTo<hasProperty,Mass>(&mahe);
    He.createRelationTo<hasProperty,Viscosity>(&muhe);
    He.createRelationTo<hasProperty,SaturationPressure>(&psathe);
    He.createRelationTo<hasProperty,SurfaceTension>(&stenhe);

    gp.createRelationTo<hasProperty,Pressure>(&p);
    gp.createRelationTo<hasProperty,Temperature>(&T);
    gp.createRelationTo<hasProperty,PressureTimeDerivative>(&dpdt);
    gp.createRelationTo<hasProperty,TemperatureTimeDerivative>(&dTdt);
    gp.createRelationTo<hasPart,HomonuclearMolecule>(&Si);
    gp.createRelationTo<hasPart,HomonuclearMolecule>(&He);
//    gp.createRelationTo<hasPart,HeteronuclearMolecule>(&CH4);

//    // test to verify if it is possible to store different entities with a vector of the class who is their ancestor
//    auto test = gp.getRelatedObject<PolyatomicEntity>();
//    for (std::size_t i = 0; i < test.size(); ++i) {
//        std::cout << test[i]->getClassName() << std::endl;
//        std::cout << test[i]->getRelatedObject<IUPAC>()[0]->data << std::endl;
//    }
//    Physical obj;

//    obj.createRelationTo<hasProperty,Pressure>(&p);
//    std::cout << obj.getRelatedObject<Pressure>()[0]->getRelatedObject<Scalar>()[0]->data;
//    std::cout << " " << obj.getRelatedObject<Pressure>()[0]->getRelatedObject<Unit>()[0]->data << std::endl;

    GasModel gm;
//    GasModelCV gm; //NOTE non rispetta la somma delle frazioni molari sempre pari a 1!!!!!!!
//    auto models = gm.isModel();
//    gm.timestep(&gp,1e-6,{-1e+30,0,0});
////    auto n = gm.get_n();

//    std::cout << "GP Temperature: " << gp.getRelatedObject<Temperature>()[0]->getRelatedObject<Scalar>()[0]->data << std::endl;
//    std::cout << "GP Pressure: " << gp.getRelatedObject<Pressure>()[0]->getRelatedObject<Scalar>()[0]->data << std::endl;
//    std::cout << "Si Molar Fraction: " << Si.getRelatedObject<MolarFraction>()[0]->getRelatedObject<Scalar>()[0]->data << std::endl;
//    std::cout << "He Molar Fraction: " << He.getRelatedObject<MolarFraction>()[0]->getRelatedObject<Scalar>()[0]->data << std::endl;
//    std::cout << "GP Average Mass: " << gm.get_average_molecular_mass(&gp) << std::endl;
//    std::cout << "Si S: " << gm.get_S(&Si,2100.) << std::endl;
//    std::cout << "Si n_sat: " << gm.get_n_sat(&Si,2100.) << std::endl;
//    std::cout <<  "GP Average Density: " << gm.get_density(&gp) << std::endl;
//    std::cout << "GP Mean Free Path: " << gm.get_mfp(&gp) << std::endl;
//    std::cout << "GP Gas Flux: " << gm.get_gas_flux(&gp) << std::endl;
//    std::cout << "GP Average Viscosity: " << gm.get_average_viscosity(&gp) << std::endl;
//    std::cout << "Si s_ten: " << Si.getRelatedObject<SurfaceTension>()[0]->get_s_ten(2000) << std::endl;

//    double Ti = gp.getRelatedObject<Temperature>()[0]->getRelatedObject<Scalar>()[0]->data;
    ClassicalNucleationTheory cnt;
//    std::cout << "Si nucleation rate: " << cnt.nucleation_rate(&Si,&gm,Ti) << std::endl;
//    std::cout << "Si stable cluster size: " << cnt.stable_cluster_size(&Si,&gm,Ti) << std::endl;
//    std::cout << "Si stable cluster diameter: " << cnt.stable_cluster_diameter(&Si,&gm,Ti) << std::endl;
//    std::cout << "Si condensation rate: " << cnt.condensation_rate(&Si,&gm,Ti) << std::endl;

    MomentModelPratsinis mm;
//    double g_cons = mm.timestep(1e-6, &gm, &cnt, &Si);
//    std::cout << "Moment model Si computed consumption: " << g_cons << std::endl;

    MomentsRun(&gm, &gp, &cnt, &mm);

    return 0;
}
