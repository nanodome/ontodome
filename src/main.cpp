#include <iostream>
#include <vector>
#include <algorithm>

#include "ontodome.h"

int main()
{
    WallClock clock;
    clock.start();

    GasMixture gas;

    MolarFraction msi(new Real(0.6), new Unit("#"));
    SingleComponentComposition si(&msi,SiliconSymbol::get_symbol());
    MolarFraction mhe(new Real(0.4), new Unit("#"));
    SingleComponentComposition he(&mhe,HeliumSymbol::get_symbol());

    Temperature T(new Real(1200.), new Unit("K"));
    Pressure p(new Real(101325.), new Unit("Pa"));
    PressureTimeDerivative dpdt(new Real(0.), new Unit("Pa/s"));
    TemperatureTimeDerivative dTdt(new Real(-1e+6), new Unit("K/s"));

    gas.createRelationsTo<hasPart,SingleComponentComposition>({&si,&he});
    gas.createRelationTo<hasProperty,Temperature>(&T);
    gas.createRelationTo<hasProperty,Pressure>(&p);
    gas.createRelationTo<hasProperty,PressureTimeDerivative>(&dpdt);
    gas.createRelationTo<hasProperty,TemperatureTimeDerivative>(&dTdt);

//	// KGDB
//	KnowledgeGeneratorsDB kgdb;
////	auto test = kgdb.ismodelfor("GasMixture");
//	auto test = kgdb.ismodelfor(new GasMixture);

//	std::cout << sisten.getLastRelation<hasModel>()->getRange()->getRelatedObjects<LatexExpression>()[0]->data << std::endl;

//	std::cout << sipsat.getLastRelation<hasModel>()->getRange()->getRelatedObjects<LatexExpression>()[0]->data << std::endl;

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
    si.getRelatedObjects<SaturationPressure>()[0]->getRelatedObjects<SaturationPressureMaterialRelation>()[0]->run();


    auto test = gas.findAll<SingleComponentComposition>();

    std::cout << "Temp: " << gas.getRelatedObjects<Temperature>()[0]->getRelatedObjects<Real>()[0]->data << std::endl;

    std::cout << "Surface Tension value: " << si.getRelatedObjects<SurfaceTension>()[0]->getRelatedObjects<Real>()[0]->data << std::endl;

    std::cout << "Saturation Pressure value: " << si.getRelatedObjects<SaturationPressure>()[0]->getRelatedObjects<Real>()[0]->data << std::endl;

    // Gas Model tests
    GasModel gm;
    gm.createRelationTo<hasModel,GasMixture>(&gas);
    double si_cons = -5e+28;
    for (int i = 1; i < 5000; i++) {
      if (i >= 1000) si_cons = 0.;
      gm.timestep(1e-7,{si_cons,0.});
    }
    gm.print();

    clock.stop();
    std::cout << "Execution time: " << clock.interval() << " s" << std::endl;

    return 0;
}
