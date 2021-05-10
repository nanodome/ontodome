#include <iostream>
#include <vector>

#include "base/thing.h"
#include "models/gasmodels/gasmodel.h"
#include "models/gasmodels/gasmodelcv.h"
#include "models/nanomodels/nucleation/cnt.h"
#include "models/nanomodels/moments/momentmodelpratsinis.h"

int main()
{
    GasMixture gp;
    Pressure p(new Scalar(101325.), new Unit("Pa"));
    Temperature T(new Scalar(500.), new Unit("K"));
    PressureTimeDerivative dpdt(new Scalar(0), new Unit("Pa/s"));
    TemperatureTimeDerivative dTdt(new Scalar(-1e+7), new Unit("K/s"));

    HomonuclearMolecule Si;
    MolarFraction msi(new Scalar(0.1), new Unit("#"));
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
    MolarFraction mch4(new Scalar(0.), new Unit("#"));
    IUPAC ch4("CH4");
    Mass mach4(new Scalar(4.002602*AMU), new Unit("kg"));
    Viscosity much4(new Scalar(5e-5), new Unit("Pa s"));
    SaturationPressure psatch4(new Vector({0,0}),new Unit("#"));
    SurfaceTension stench4(new Vector({0.,0.,0.}),new Unit("#"));
    CH4.createRelationTo<hasProperty,MolarFraction>(&mch4);
    CH4.createRelationTo<hasProperPart,IUPAC>(&ch4);
    CH4.createRelationTo<hasProperty,Mass>(&mach4);
    CH4.createRelationTo<hasProperty,Viscosity>(&much4);
    CH4.createRelationTo<hasProperty,SaturationPressure>(&psatch4);
    CH4.createRelationTo<hasProperty,SurfaceTension>(&stench4);

    HomonuclearMolecule He;
    MolarFraction mhe(new Scalar(0.9), new Unit("#"));
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
    gp.createRelationTo<hasPart,HeteronuclearMolecule>(&CH4);

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
    auto models = gm.isModel();
    gm.run(&gp,1e-6,{-1e+30,0,0});
//    auto n = gm.get_n();

    std::cout << "GP Temperature: " << gp.getRelatedObject<Temperature>()[0]->getRelatedObject<Scalar>()[0]->data << std::endl;
    std::cout << "GP Pressure: " << gp.getRelatedObject<Pressure>()[0]->getRelatedObject<Scalar>()[0]->data << std::endl;
    std::cout << "Si Molar Fraction: " << Si.getRelatedObject<MolarFraction>()[0]->getRelatedObject<Scalar>()[0]->data << std::endl;
    std::cout << "He Molar Fraction: " << He.getRelatedObject<MolarFraction>()[0]->getRelatedObject<Scalar>()[0]->data << std::endl;
    std::cout << "GP Average Mass: " << gm.get_average_molecular_mass(&gp) << std::endl;
    std::cout << "Si S: " << gm.get_S(&Si,2100.) << std::endl;
    std::cout << "Si n_sat: " << gm.get_n_sat(&Si,2100.) << std::endl;
    std::cout <<  "GP Average Density: " << gm.get_density(&gp) << std::endl;
    std::cout << "GP Mean Free Path: " << gm.get_mfp(&gp) << std::endl;
    std::cout << "GP Gas Flux: " << gm.get_gas_flux(&gp) << std::endl;
    std::cout << "GP Average Viscosity: " << gm.get_average_viscosity(&gp) << std::endl;
    std::cout << "Si s_ten: " << Si.getRelatedObject<SurfaceTension>()[0]->get_s_ten(2000) << std::endl;

    double Ti = gp.getRelatedObject<Temperature>()[0]->getRelatedObject<Scalar>()[0]->data;
    ClassicalNucleationTheory cnt;
    std::cout << "Si nucleation rate: " << cnt.nucleation_rate(&Si,&gm,Ti) << std::endl;
    std::cout << "Si stable cluster size: " << cnt.stable_cluster_size(&Si,&gm,Ti) << std::endl;
    std::cout << "Si stable cluster diameter: " << cnt.stable_cluster_diameter(&Si,&gm,Ti) << std::endl;
    std::cout << "Si condensation rate: " << cnt.condensation_rate(&Si,&gm,Ti) << std::endl;

    MomentModelPratsinis mm;
    double g_cons = mm.timestep(1e-6, &gm, &cnt, &Si);
    std::cout << "Moment model Si computed consumption: " << g_cons << std::endl;


    double te = 0;
}
