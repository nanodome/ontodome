#include <iostream>
#include <vector>

#include "base/thing.h"
#include "models/gasmodel.h"
#include "models/gasmodelcv.h"

int main()
{
    GasMixture gp;
    Pressure p(new Scalar(101325.), new Unit("Pa"));
    Temperature T(new Scalar(300.), new Unit("K"));
    PressureTimeDerivative dpdt(new Scalar(-1e+7), new Unit("Pa/s"));
    TemperatureTimeDerivative dTdt(new Scalar(-1e+7), new Unit("K/s"));

    HomonuclearMolecule Si;
    MolarFraction msi(new Scalar(0.1), new Unit("#"));
    IUPAC si("Si");
    Mass masi(new Scalar(28.085*AMU), new Unit("kg"));
    Viscosity musi(new Scalar(7e-5), new Unit("Pa s"));
    SaturationPressure psatsi(new Vector({7.5341,23399.}),new Unit("#"));
    SurfaceTension stensi(new Vector({0.732,0.000086,1685.}),new Unit("#"));
    Si.createRelationTo<hasProperty,MolarFraction>(&msi);
    Si.createRelationTo<hasProperty,IUPAC>(&si);
    Si.createRelationTo<hasProperty,Mass>(&masi);
    Si.createRelationTo<hasProperty,Viscosity>(&musi);
    Si.createRelationTo<hasProperty,SaturationPressure>(&psatsi);
    Si.createRelationTo<hasProperty,SurfaceTension>(&stensi);

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

//    Physical obj;

//    obj.createRelationTo<hasProperty,Pressure>(&p);
//    std::cout << obj.getRelatedObject<Pressure>()[0]->getRelatedObject<Scalar>()[0]->data;
//    std::cout << " " << obj.getRelatedObject<Pressure>()[0]->getRelatedObject<Unit>()[0]->data << std::endl;

//    GasModel gm;
    GasModelCV gm; //NOTE non rispetta la somma delle frazioni molari sempre pari a 1!!!!!!!
    auto models = gm.isModel();
    gm.run(&gp,1e-6,{-1e+30,0});
    auto n = gm.get_n();

    std::cout << gp.getRelatedObject<Temperature>()[0]->getRelatedObject<Scalar>()[0]->data << std::endl;
    std::cout << gp.getRelatedObject<Pressure>()[0]->getRelatedObject<Scalar>()[0]->data << std::endl;
    std::cout << Si.getRelatedObject<MolarFraction>()[0]->getRelatedObject<Scalar>()[0]->data << std::endl;
    std::cout << He.getRelatedObject<MolarFraction>()[0]->getRelatedObject<Scalar>()[0]->data << std::endl;
    std::cout << gm.get_average_molecular_mass(&gp) << std::endl;
    std::cout << gm.get_S<HomonuclearMolecule>(&Si,2100.) << std::endl;
    std::cout << gm.get_n_sat<HomonuclearMolecule>(&Si,2100.) << std::endl;
    std::cout << gm.get_density(&gp) << std::endl;
    std::cout << gm.get_average_molecular_mass(&gp) << std::endl;
    std::cout << gm.get_mfp(&gp) << std::endl;
    std::cout << gm.get_gas_flux(&gp) << std::endl;
    std::cout << gm.get_average_viscosity(&gp) << std::endl;

    double te = 0;
}
