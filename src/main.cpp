#include <iostream>
#include <vector>
#include <algorithm>

#include "ontodome.h"

//#include "base/thing.h"
//#include "knowledgegenerators/surfacetensionmodels/surfacetensionpolynomialmodel.h"

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

    gas.createRelationsTo<hasPart,SingleComponentComposition>({&si,&he});
    gas.createRelationTo<hasProperty,Temperature>(&T);

//	// KGDB
//	KnowledgeGeneratorsDB kgdb;
////	auto test = kgdb.ismodelfor("GasMixture");
//	auto test = kgdb.ismodelfor(new GasMixture);

//	std::cout << sisten.getLastRelation<hasModel>()->getRange()->getRelatedObjects<LatexExpression>()[0]->data << std::endl;

//	std::cout << sipsat.getLastRelation<hasModel>()->getRange()->getRelatedObjects<LatexExpression>()[0]->data << std::endl;

    SurfaceTensionPolynomialModel stpm;
    SurfaceTensionMaterialRelation stmr;
    SurfaceTension st(new Real(0.), new Unit("N/m"));

    stmr.createRelationTo<hasSoftwareModel,SurfaceTensionPolynomialModel>(&stpm);
    st.createRelationTo<hasMathematicalModel,SurfaceTensionMaterialRelation>(&stmr);
    si.createRelationTo<hasProperty,SurfaceTension>(&st);
    si.getRelatedObjects<SurfaceTension>()[0]->getRelatedObjects<MaterialRelation>()[0]->run();

    auto test = gas.findAll<SingleComponentComposition>();

    std::cout << "Temp: " << gas.getRelatedObjects<Temperature>()[0]->getRelatedObjects<Real>()[0]->data << std::endl;

    std::cout << "Surface Tension value: " << si.getRelatedObjects<SurfaceTension>()[0]->getRelatedObjects<Real>()[0]->data << std::endl;

    clock.stop();
    std::cout << "Execution time: " << clock.interval() << " s" << std::endl;

    return 0;
}
