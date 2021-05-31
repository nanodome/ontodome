#ifndef SURFACETENSIONMATERIALRELATION_H
#define SURFACETENSIONMATERIALRELATION_H

//#include "knowledgegenerators.h"
#include "../base/thing.h"

class SurfaceTensionMaterialRelation : public MaterialRelation {

public:
    std::string getClassName() const { return "SurfaceTensionMaterialRelation"; }

    void run()
    {
        // Find and run first related Software object model
        // std::cout << "Software Model is: " << find<SoftwareModel>()->getClassName() << std::endl;
        find<SoftwareModel>()->run();
    }
};

#endif // SURFACETENSIONMATERIALRELATION_H
