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

    return 0;
}
