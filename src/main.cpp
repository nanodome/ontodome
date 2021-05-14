#include <iostream>
#include <vector>

#include "ontodome.h"
#include "base/thing.h"
/*#include "models/gasmodels/gasmodel.h"
#include "models/gasmodels/gasmodelcv.h"
#include "models/nanomodels/nucleation/cnt.h"
#include "models/nanomodels/moments/momentmodelpratsinis.h"
#include "models/statemodels/stateinterpolator.h"

/*void MomentsRun(GasModel* gm, GasMixture* GP, ClassicalNucleationTheory* cnt, MomentModelPratsinis* mm) {

  // Condensing species datas
  // Get the initial state of the domain
  auto gp = GP->get_State(0.);

  // Get the initial parameters
  auto spec = gp->getRelatedObject<PolyatomicEntity>()[0];
  std::valarray<double> w_cons;
  w_cons.resize(gp->getRelatedObject<PolyatomicEntity>().size());
  auto spec_name = spec->getRelatedObject<IUPAC>()[0]->data;

  // Initialize the models
  gm->initialize(gp);
  cnt->initialize(spec);
  mm->initialize(cnt,spec);

  // simulation settings
//  const double dt = 1e-8;
  const double dt = 5e-7;
  double t = 0.;
  const double t_end = 1.;
  int iter = 0;

  const int PRINT_EVERY = 100;
  const double SAVE_STATE_EVERY = 5e-3;
  double last_save = 0;

  // loop over timesteps until the final temperature
  while(t <= t_end) {

      double T = gm->get_T();
      double p = gm->get_p();
      double S = gm->get_S(spec,T);
//      double J = cnt->nucleation_rate(spec,gm,T);
//      double j = cnt->stable_cluster_size(spec,gm,T);

      if (T < 300. || p < 10.) { break; }
//      if (T <= 300.) {
//          gm->dTdt = 0;
//      }

      // moment method timestep and species consumption retrieval
      w_cons[0] = - mm->timestep(dt, gm, cnt, spec);

      // updating the gas phase
      gm->timestep(dt,w_cons);

      ++iter;
      t+=dt;

      if (counter_trigger(iter, PRINT_EVERY)) {
          gm->print();

          std::cout
              << "Time= "<< t << '\t'
              << "Temp= " << T << '\t'
              << "Sat_" << spec_name << "= " << S << '\t'
              << "w_" << spec_name << "= " << gm->w[0] << '\t'
              << spec_name << "_cons= " << -w_cons[0] << '\t'
              << "Mean Diam_" << spec_name << "= " << mm->get_mean_diameter() << '\t'
              << "M0_" << spec_name << "= " << mm->get_M0() << '\t'
              << "M1_" << spec_name << "= " << mm->get_M1() << '\t'
              << "M2_" << spec_name << "= " << mm->get_M2() << '\t'
              << std::endl << std::endl;
      }

      // Save the GasMixture state every SAVE_STATE_EVERY [s]
      if ((t-last_save) > SAVE_STATE_EVERY) {
        gm->add_temporal_state(GP,t);
        last_save = t;
      }
  }

  // add the final state of GasMixutre
  gm->add_temporal_state(GP,t);
}*/

int main()
{
    WallClock clock;
    clock.start();

    double ww = 0.0005;
    double T_start = 4000.;

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

    HeteronuclearMolecule CH4;
    MolarFraction mch4(new Scalar((1-ww)*0.025), new Unit("#"));
    IUPAC ch4("CH4");
    Mass mach4(new Scalar(16.04*AMU), new Unit("kg"));
    Viscosity much4(new Scalar(5e-5), new Unit("Pa s"));
    SaturationPressure psatch4(new Vector({0,0}),new Unit("Pa"));
    SurfaceTension stench4(new Vector({0.,0.,0.}),new Unit("N/m"));
    CH4.createRelationTo<hasProperty,MolarFraction>(&mch4);
    CH4.createRelationTo<hasProperPart,IUPAC>(&ch4);
    CH4.createRelationTo<hasProperty,Mass>(&mach4);
    CH4.createRelationTo<hasProperty,Viscosity>(&much4);
    CH4.createRelationTo<hasProperty,SaturationPressure>(&psatch4);
    CH4.createRelationTo<hasProperty,SurfaceTension>(&stench4);

    HomonuclearMolecule H;
    MolarFraction mh(new Scalar((1-ww)*0.025), new Unit("#"));
    IUPAC h("H");
    Mass mah(new Scalar(1.0006*AMU), new Unit("kg"));
    Viscosity muh(new Scalar(5e-5), new Unit("Pa s"));
    SaturationPressure psath(new Vector({0,0}),new Unit("Pa"));
    SurfaceTension stenh(new Vector({0.,0.,0.}),new Unit("N/m"));
    H.createRelationTo<hasProperty,MolarFraction>(&mh);
    H.createRelationTo<hasProperPart,IUPAC>(&h);
    H.createRelationTo<hasProperty,Mass>(&mah);
    H.createRelationTo<hasProperty,Viscosity>(&muh);
    H.createRelationTo<hasProperty,SaturationPressure>(&psath);
    H.createRelationTo<hasProperty,SurfaceTension>(&stenh);

    HomonuclearMolecule He;
    MolarFraction mhe(new Scalar((1-ww)*0.95), new Unit("#"));
    IUPAC he("He");
    Mass mahe(new Scalar(4.002602*AMU), new Unit("kg"));
    Viscosity muhe(new Scalar(5e-5), new Unit("Pa s"));
    SaturationPressure psathe(new Vector({0,0}),new Unit("Pa"));
    SurfaceTension stenhe(new Vector({0.,0.,0.}),new Unit("N/m"));
    He.createRelationTo<hasProperty,MolarFraction>(&mhe);
    He.createRelationTo<hasProperty,IUPAC>(&he);
    He.createRelationTo<hasProperty,Mass>(&mahe);
    He.createRelationTo<hasProperty,Viscosity>(&muhe);
    He.createRelationTo<hasProperty,SaturationPressure>(&psathe);
    He.createRelationTo<hasProperty,SurfaceTension>(&stenhe);

 /*   GasMixture gp;
    Pressure p(new Scalar(101325.), new Unit("Pa"));
    Temperature T(new Scalar(T_start), new Unit("K"));
    PressureTimeDerivative dpdt(new Scalar(0.), new Unit("Pa/s"));
    TemperatureTimeDerivative dTdt(new Scalar(-1e+5), new Unit("K/s"));

    gp.initialize<hasProperty,Quantity>({&p,&T,&dpdt,&dTdt});
    gp.initialize<hasPart,PolyatomicEntity>({&Si,&He,&CH4,&H});

    GasModel gm;
//    GasModelCV gm;
//    auto models = gm.isModel();

    ClassicalNucleationTheory cnt;

    MomentModelPratsinis mm;

    MomentsRun(&gm, &gp, &cnt, &mm);

    StateInterpolator interp;

    auto test = interp.intepolate_state(&gp,0.0124);
    std::cout << "Interpolated temperature is: " << test->getRelatedObject<Temperature>()[0]->getRelatedObject<Scalar>()[0]->data << std::endl;
    std::cout << "Interpolated pressure is: " << test->getRelatedObject<Pressure>()[0]->getRelatedObject<Scalar>()[0]->data << std::endl;

    clock.stop();
    std::cout << "Execution time: " << clock.interval() << " s" << std::endl;

    return 0;*/
}
