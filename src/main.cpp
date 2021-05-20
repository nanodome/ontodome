#include <iostream>
#include <vector>
#include <algorithm>

#include "ontodome.h"

//void MomentsRun(GasModel* gm, GasMixture* GP, ClassicalNucleationTheory* cnt, MomentModelPratsinis* mm) {

//  // Condensing species datas
//  // Get the initial state of the domain
//  auto gp = GP->get_State(0.);

//  // Get the initial parameters
//  auto spec = gp->getRelatedObject<PolyatomicEntity>()[0];
//  std::valarray<double> w_cons;
//  w_cons.resize(gp->getRelatedObject<PolyatomicEntity>().size());
//  auto spec_name = spec->getRelatedObject<IUPAC>()[0]->data;

//  // Initialize the models
//  gm->initialize(gp);
//  cnt->initialize(spec);
//  mm->initialize(cnt,spec);

//  // simulation settings
////  const double dt = 1e-8;
//  const double dt = 5e-7;
//  double t = 0.;
//  const double t_end = 1.;
//  int iter = 0;

//  const int PRINT_EVERY = 100;
//  const double SAVE_STATE_EVERY = 5e-3;
//  double last_save = 0;

//  // loop over timesteps until the final temperature
//  while(t <= t_end) {

//      double T = gm->get_T();
//      double p = gm->get_p();
//      double S = gm->get_S(spec,T);
////      double J = cnt->nucleation_rate(spec,gm,T);
////      double j = cnt->stable_cluster_size(spec,gm,T);

//      if (T < 300. || p < 10.) { break; }
////      if (T <= 300.) {
////          gm->dTdt = 0;
////      }

//      // moment method timestep and species consumption retrieval
//      w_cons[0] = - mm->timestep(dt, gm, cnt, spec);

//      // updating the gas phase
//      gm->timestep(dt,w_cons);

//      ++iter;
//      t+=dt;

//      if (counter_trigger(iter, PRINT_EVERY)) {
//          gm->print();

//          std::cout
//              << "Time= "<< t << '\t'
//              << "Temp= " << T << '\t'
//              << "Sat_" << spec_name << "= " << S << '\t'
//              << "w_" << spec_name << "= " << gm->w[0] << '\t'
//              << spec_name << "_cons= " << -w_cons[0] << '\t'
//              << "Mean Diam_" << spec_name << "= " << mm->get_mean_diameter() << '\t'
//              << "M0_" << spec_name << "= " << mm->get_M0() << '\t'
//              << "M1_" << spec_name << "= " << mm->get_M1() << '\t'
//              << "M2_" << spec_name << "= " << mm->get_M2() << '\t'
//              << std::endl << std::endl;
//      }

//      // Save the GasMixture state every SAVE_STATE_EVERY [s]
//      if ((t-last_save) > SAVE_STATE_EVERY) {
//        gm->add_temporal_state(GP,t);
//        last_save = t;
//      }
//  }

//  // add the final state of GasMixutre
//  gm->add_temporal_state(GP,t);
//}

int main()
{
	WallClock clock;
	clock.start();

//	double ww = 0.0005;

//        HomonuclearMolecule Si;
//        MolarFraction msi(new Scalar(ww), new Unit("#"));
//        IUPAC si("Si");
//        Mass masi(new Scalar(28.085*AMU), new Unit("kg"));
//        Viscosity musi(new Scalar(7e-5), new Unit("Pa s"));
//        BulkDensityLiquid bdlsi(new Scalar(2570.), new Unit("kg/m3"));
//        BulkDensitySolid bdssi(new Scalar(2329.), new Unit("kg/m3"));
//        MeltingPoint mpsi(new Scalar(1687.), new Unit("K"));
//        SaturationPressure psatsi(new Scalar(0.), new Unit("Pa"));
//        SurfaceTension stensi(new Scalar(0.), new Unit("N/m"));

//        Si.createRelationsTo<hasProperty,Thing>({&msi,&si,&masi,&musi,&bdlsi,&bdssi,&mpsi,&psatsi,&stensi});

////        SaturationPressureModel sipsat({7.5341,23399.});
//        SaturationPressureModel sipsat(Si.getLastRelation<IUPAC>());
//        psatsi.createRelationTo<hasModel,SaturationPressureModel>(&sipsat);

//        SurfaceTensionModel sisten({0.732,0.000086,1685.});
//        stensi.createRelationTo<hasModel,SurfaceTensionModel>(&sisten);

//        Gas gasgas;
//        SingleComponentComposition aa(&mch4,Methane::get_symbol());
//        MolarFraction mch4(new Scalar((1-ww)*0.025), new Unit("#"));
//        SingleComponentComposition aaa(&mhe,Helium::get_symbol());
//        MolarFraction mch4(new Scalar((1-ww)*0.025), new Unit("#"));
//        IUPAC ch4("CH4");
//        Mass mach4(new Scalar(16.04*AMU), new Unit("kg"));
//        Viscosity much4(new Scalar(5e-5), new Unit("Pa s"));
//        SaturationPressure psatch4(new Scalar(0),new Unit("Pa"));
//        SurfaceTension stench4(new Scalar(0),new Unit("N/m"));

//        CH4.createRelationsTo<hasProperty,Thing>({&mch4,&ch4,&mach4,&much4,&psatch4,&stench4});

//        GasMixture gp;

//	GasMixture gp0;
//	Pressure p0(new Scalar(101325.), new Unit("Pa"));
//	Temperature T0(new Scalar(1700), new Unit("K"));
//	PressureTimeDerivative dpdt0(new Scalar(0.), new Unit("Pa/s"));
//	TemperatureTimeDerivative dTdt0(new Scalar(-1e+5), new Unit("K/s"));
//	Time t0(new Scalar(0.), new Unit("s"));

//	gp0.createRelationsTo<hasProperty,Thing>({&p0,&T0,&dpdt0,&dTdt0,&t0});
//	gp0.createRelationsTo<hasPart,Thing>({&Si,&He,&CH4,&H});
//	gp.createRelationTo<hasTemporalDirectPart,Thing>(&gp0);

//	// SaturationPressure and SurfaceTension models tests
//	double TT = gp.getLastRelation<hasTemporalDirectPart>()->getRange()->getRelatedScalarObjects<Temperature>()[0];

//	auto sat = Si1.getRelatedObjects<SaturationPressure>()[0];
//	sat->getRelatedObjects<SaturationPressureModel>()[0]->run(TT);
//	std::cout << "SaturationPressure at T = " << TT << " is: " << Si1.getRelatedScalarObjects<SaturationPressure>()[0] << std::endl;

//	auto sten = Si1.getRelatedObjects<SurfaceTension>()[0];
//	sten->getRelatedObjects<SurfaceTensionModel>()[0]->run(TT);
//	std::cout << "SurfaceTension at T = " << TT << " is: " << Si1.getRelatedScalarObjects<SurfaceTension>()[0] << std::endl;

//	// KGDB
//	KnowledgeGeneratorsDB kgdb;
////	auto test = kgdb.ismodelfor("GasMixture");
//	auto test = kgdb.ismodelfor(new GasMixture);

//	std::cout << sisten.getLastRelation<hasModel>()->getRange()->getRelatedObjects<LatexExpression>()[0]->data << std::endl;

//	std::cout << sipsat.getLastRelation<hasModel>()->getRange()->getRelatedObjects<LatexExpression>()[0]->data << std::endl;

	clock.stop();
	std::cout << "Execution time: " << clock.interval() << " s" << std::endl;

	return 0;
}
