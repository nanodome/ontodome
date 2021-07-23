#include "ontodome.h"

int main()
{
    // Setting a clock to keep track of computational time
    WallClock clock;
    clock.start();

    // GasMixture initial state declaration
    GasMixture gas;

    Time t(new Real(0), new Unit("s"));
    gas.createRelationTo<hasProperty,Time>(&t);

    MolarFraction msi(new Real(0.05), new Unit("#"));
    SingleComponentComposition si(&msi,SiliconSymbol::get_symbol());
    MolarFraction mhe(new Real(0.95), new Unit("#"));
    SingleComponentComposition he(&mhe,HeliumSymbol::get_symbol());

    Temperature T(new Real(4000.), new Unit("K"));
    Pressure p(new Real(101325.), new Unit("Pa"));
    PressureTimeDerivative dpdt(new Real(0.), new Unit("Pa/s"));
    TemperatureTimeDerivative dTdt(new Real(-1e+6), new Unit("K/s"));

    gas.createRelationsTo<hasPart,SingleComponentComposition>({&si,&he});
    gas.createRelationTo<hasProperty,Temperature>(&T);
    gas.createRelationTo<hasProperty,Pressure>(&p);
    gas.createRelationTo<hasProperty,PressureTimeDerivative>(&dpdt);
    gas.createRelationTo<hasProperty,TemperatureTimeDerivative>(&dTdt);

    // Surface Tension and Saturation Pressure settings
    SurfaceTensionPolynomialSoftwareModel stpm;
    SurfaceTensionMaterialRelation stmr;
    SurfaceTension st(new Real(0.), new Unit("N/m"));

    stmr.createRelationTo<hasSoftwareModel,SurfaceTensionPolynomialSoftwareModel>(&stpm);
    st.createRelationTo<hasMathematicalModel,SurfaceTensionMaterialRelation>(&stmr);
    si.createRelationTo<hasProperty,SurfaceTension>(&st);
    si.getRelatedObjects<SurfaceTension>()[0]->getRelatedObjects<SurfaceTensionMaterialRelation>()[0]->run();

    SaturationPressurePolynomialSoftwareModel sapm;
    SaturationPressureMaterialRelation samr;
    SaturationPressure sa(new Real(0.), new Unit("Pa"));

    samr.createRelationTo<hasSoftwareModel,SaturationPressurePolynomialSoftwareModel>(&sapm);
    sa.createRelationTo<hasMathematicalModel,SaturationPressureMaterialRelation>(&samr);
    si.createRelationTo<hasProperty,SaturationPressure>(&sa);

    // Si properties
    Mass sim(new Real(28.0855*AMU), new Unit("kg"));
    sim.createRelationTo<hasScalarProperty,SingleComponentComposition>(&si);

    BulkDensityLiquid sibl(new Real(2570.), new Unit("kg/m3"));
    sibl.createRelationTo<hasScalarProperty,SingleComponentComposition>(&si);

    BulkDensitySolid sibs(new Real(2329.), new Unit("kg/m3"));
    sibs.createRelationTo<hasScalarProperty,SingleComponentComposition>(&si);

    MeltingPoint melp(new Real(1687.), new Unit("K"));
    melp.createRelationTo<hasScalarProperty,SingleComponentComposition>(&si);

    // Models settings
    Time dt(new Real(1e-7), new Unit("s")); // simulation timestep. This definition is not ontologically correct

    GasModel gm;
    gm.createRelationTo<hasModel,GasMixture>(&gas);

    ClassicalNucleationTheory cnt;
    cnt.createRelationTo<hasModel,SingleComponentComposition>(&si);

    MomentModelPratsinis mom;
    mom.createRelationTo<hasModel,SingleComponentComposition>(&si);

    PBMFractalParticlePhase<PBMAggregate<Particle>> pp(1.61, 5.0e-16);
    pp.createRelationTo<hasModel,SingleComponentComposition>(&si);
    // set max and min number of aggregates
    pp.set_max_aggregates(2000);
    pp.set_min_aggregates(1990);

    // Example of usage for the moments method
    int PRINT_EVERY = 1000;
    int iter = 0;

    // PBM loop
    // loop over timesteps
    while(*t.onData() < 0.002) {

      if (gm.get_T() < 300)
      { *dTdt.onData() = 0; }

      // species source term for the gas phase
      double g_si = 0.0;

      // calculate the timestep using an exponential waiting time
      double R_tot = pp.get_total_processes_rate(&gm,&cnt);
      double rho = ndm::uniform_double_distr(ndm::rand_gen);

      // exponential waiting time
      double dt = -log(rho)/R_tot;

      // Strang first step
      gm.timestep(dt/2.0,{0.0,0.0});
      pp.volume_expansion(dt/2.0,&gm);

      // Strang second step
      g_si += pp.timestep(dt,&gm,&cnt,&si);
      gm.timestep(dt,{-g_si,0.0});

      // Strang third step
      gm.timestep(dt/2.0,{0.0,0.0});
      pp.volume_expansion(dt/2.0,&gm);

  //            *t.onData() += *dt.onData();
      iter++;
      if(counter_trigger(iter,PRINT_EVERY)) {

          clock.stop();

          std::cout << *t.onData() << '\t'				// time
                    << gm.get_T() << '\t'					// temperature
                    << gm.get_S(&si) << '\t'				// supersaturation (S)
                    << cnt.nucleation_rate() << '\t'			// J
                    << gm.get_n() << '\t'					// ns
                    << cnt.stable_cluster_diameter() << '\t'		// j
                    << pp.get_mean_particles_number() << '\t'		// N_m
                    << pp.get_mean_sintering_level() << '\t'		//
                    << pp.get_aggregates_mean_spherical_diameter() << '\t'
                    << pp.get_aggregates_number() << '\t'
                    << pp.get_aggregates_density() << '\t'
                    << pp.get_volume() << '\t'
                    << pp.get_mean_fractal_dimension() << '\t'
                    << clock.interval()/PRINT_EVERY << std::endl;

          clock.start();
      }
      }


/*    while ( *t.onData() <= 0.02) {
      if ( gm.get_T() < 600.) {
        *dTdt.onData() = 0.;
      }

      double g_cons = mom.timestep(*dt.onData());

      gm.timestep(*dt.onData(), {g_cons, 0.0});

      *t.onData() += *dt.onData();
      iter += 1;

      if (counter_trigger(iter, PRINT_EVERY)) {

          auto spec_name = si.name;
          std::cout
              << "Time= "<< *t.onData() << '\t'
              << "Temp= " << gm.get_T() << '\t'
              << "Sat_" << spec_name << "= " << gm.get_S(&si) << '\t'
              << spec_name << "_cons= " << g_cons << '\t'
              << "Mean Diam_" << spec_name << "= " << mom.get_mean_diameter() << '\t'
              << "M0_" << spec_name << "= " << mom.get_n_density() << '\t'
              << "M1_" << spec_name << "= " << mom.get_M1() << '\t'
              << "M2_" << spec_name << "= " << mom.get_M2() << '\t'
              << std::endl << std::endl;
      }
    }
*/
    gm.print();

    clock.stop();
    std::cout << "Execution time: " << clock.interval() << " s" << std::endl;

    return 0;
}
