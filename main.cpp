#include "ontodome.h"

//void particle() {

//	/*	On a multivariate population balance model to describe the structure and composition of silica nanoparticles.
//		Shraddha Shekara, William J. Menza, Alastair J. Smith, Markus Kraft, Wolfgang Wagnerb
//		Computers and Chemical Engineering	*/

//    WallClock clock;

//    Species si("Si");
//    Species ar("Ar");

//    // create the vector of species for the gas phase
//    std::vector<Species> gas_species = {si, ar};

//    // gas phase info
//    double p = 1.01e5; // pressure [Pa]
//    double dTdt = -1.0e7; // temperature gradient [K/s]

//    double T_start = 3000.0;

//    // set pressure, temperature, species and initial molar concentration
//    GasPhaseCV gp(p, T_start, gas_species, {0.01, 0.99});

//    // setup the nucleation theory we want to use
//    ClassicalNucleationTheory cnt(si);

//    // setup the particle phase we want to use
//    PBMFractalParticlePhase<PBMAggregate<Particle>> pp(1.61, 5.0e-16);

//	// set max and min number of aggregates
//	pp.set_max_aggregates(2000);
//	pp.set_min_aggregates(1990);

//    double t = 0.0;
//    int iter = 0;

//    const int PRINT_EVERY = 1000;
//	const int PRINT_HEADLINE = PRINT_EVERY * 5;


//	const int SAVE_EVERY = 1000;
//	std::ofstream plot_data;
//	std::string filename = "PBM_plot.dat";
//	plot_data.open(filename);
//	// Print headlines
//	plot_data << "Time[sec]" << '\t'
//		<< "Temp[K]" << '\t'
//		<< "SuperSaturation ratio" << '\t'
//		<< "Nucl Rate" << '\t'
//		<< "Species # density" << '\t'
//		<< "Stable cluster size[m]" << '\t'
//		<< "AVG Part Num[#]: " << '\t'
//		<< "Sint level[%]" << '\t'
//		<< "AVG diameter[m]" << '\t'
//		<< "Agg. #[#]" << '\t'
//		<< "Agg density[#/m3]" << '\t'
//		<< "Volume[m]" << '\t'
//		<< "AVG fract dim" << '\t'
//		<< "ts exec time"
//		<< std::endl;

//	// PSD Data files
//	const double PSD_DATA = 1.0e-5;
//	double save_step = 0.0;
//	// Prinary particles
//	std::ofstream part_sizes_file;
//	std::string part_size_filename = "particles_sizes.dat";
//	part_sizes_file.open(part_size_filename);
//	// Aggregates diameters
//	std::ofstream agg_sizes_file;
//	std::string agg_size_filename = "aggregates_sizes.dat";
//	agg_sizes_file.open(agg_size_filename);

//    clock.start();

//    // loop over timesteps
//    while(t < 0.002) {

//		if (gp.get_T() < 300)
//		{ dTdt = 0; }

//        // species source term for the gas phase
//        double g_si = 0.0;

//        // calculate the timestep using an exponential waiting time
//        double R_tot = pp.get_total_processes_rate(gp,cnt);
//        double rho = ndm::uniform_double_distr(ndm::rand_gen);

//        // exponential waiting time
//        double dt = -log(rho)/R_tot;

//        // Strang first step
//        gp.timestep(dt/2.0,dTdt,0);
//        pp.volume_expansion(dt/2.0,gp);

//        // Strang second step
//        g_si += pp.timestep(dt,gp,cnt);
//        gp.timestep(dt,0,0,{-g_si,0.0});

//        // Strang third step
//        gp.timestep(dt/2.0,dTdt,0);
//        pp.volume_expansion(dt/2.0,gp);

//        t += dt;
//		save_step += dt;
//        iter++;

//        if(counter_trigger(iter,PRINT_EVERY)) {

//            clock.stop();

//            double T = gp.get_T();
//            double S = gp.get_S("Si");

//            std::cout << t << '\t'									// time
//                      << gp.get_T() << '\t'							// temperature
//                      << gp.get_S("Si") << '\t'						// supersaturation (S)
//                      << cnt.nucleation_rate(T, S) << '\t'			// J
//                      << gp.get_n("Si") << '\t'						// ns
//                      << cnt.stable_cluster_size(T, S) << '\t'		// j
//                      << pp.get_mean_particles_number() << '\t'		// N_m
//                      << pp.get_mean_sintering_level() << '\t'		//
//                      << pp.get_aggregates_mean_spherical_diameter() << '\t'
//                      << pp.get_aggregates_number() << '\t'
//                      << pp.get_aggregates_density() << '\t'
//					  << pp.get_volume() << '\t'
//					  << pp.get_mean_fractal_dimension() << '\t'
//                      << clock.interval()/PRINT_EVERY << std::endl;

//            clock.start();
//        }

//		// Print file for plotting
//		if (counter_trigger(iter, SAVE_EVERY)) {

//			clock.stop();

//			double T = gp.get_T();
//			double S = gp.get_S("Si");

//			plot_data
//				<< t << '\t'
//				<< gp.get_T() << '\t'
//				<< gp.get_S("Si") << '\t'
//				<< cnt.nucleation_rate(T, S) << '\t'
//				<< gp.get_n("Si") << '\t'
//				<< cnt.stable_cluster_size(T, S) << '\t'
//				<< pp.get_mean_particles_number() << '\t'
//				<< pp.get_mean_sintering_level() << '\t'
//				<< pp.get_aggregates_mean_spherical_diameter() << '\t'
//				<< pp.get_aggregates_number() << '\t'
//				<< pp.get_aggregates_density() << '\t'
//				<< pp.get_volume() << '\t'
//				<< pp.get_mean_fractal_dimension() << '\t'
//				<< clock.interval() / PRINT_EVERY << std::endl;

//			clock.start();
//		}

//		// Save particles sizes for PSD
//		//if (save_step>=PSD_DATA) {

//		//	clock.stop();

//		//	std::valarray<double> particles_sizes = pp.get_particles_sizes();
//		//	std::valarray<double> aggregates_sizes = pp.get_aggregates_sizes();

//		//	// print particles sizes
//		//	for (int i = 0; i < particles_sizes.size(); i++) {
//		//		part_sizes_file << particles_sizes[i] << " ";
//		//	}
//		//	part_sizes_file << std::endl;
//		//	// print aggregates sizes
//		//	for (int i = 0; i < aggregates_sizes.size(); i++) {
//		//		agg_sizes_file << aggregates_sizes[i] << " ";
//		//	}
//		//	agg_sizes_file << std::endl;

//		//	save_step = 0.0;

//		//	clock.start();
//		//}
//    }
//	// Close plot file
//	plot_data.close();

//	// Anyway at the end save aggregates data for PSD
//	std::valarray<double> particles_sizes = pp.get_particles_sizes();
//	std::valarray<double> aggregates_sizes = pp.get_aggregates_sizes();

//	// print particles sizes
//	for (int i = 0; i < particles_sizes.size(); i++) {
//		part_sizes_file << particles_sizes[i] << " ";
//	}
//	part_sizes_file << std::endl;
//	// print aggregates sizes
//	for (int i = 0; i < aggregates_sizes.size(); i++) {
//		agg_sizes_file << aggregates_sizes[i] << " ";
//	}
//	agg_sizes_file << std::endl;

//	part_sizes_file.close();
//	agg_sizes_file.close();

//}

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

	// Example of usage for the moments method
	double* ts = t.onData();
	int PRINT_EVERY = 1000;
	int iter = 0;

	while ( *ts <= 0.02) {
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
			  << "Time= "<< *ts << '\t'
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

	gm.print();

	clock.stop();
	std::cout << "Execution time: " << clock.interval() << " s" << std::endl;

	gas.findNearest<GasModel>();

	return 0;
}
