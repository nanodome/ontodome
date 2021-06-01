#ifndef SATURATIONPRESSUREMATERIALRELATION_H
#define SATURATIONPRESSUREMATERIALRELATION_H

//#include "knowledgegenerators.h"
#include "../base/thing.h"

class SaturationPressureMaterialRelation : public MaterialRelation {

public:
    std::string getClassName() const { return "SaturationPressureMaterialRelation"; }

    void run()
    {
        // Find and run first related Software object model
        // std::cout << "Software Model is: " << find<SoftwareModel>()->getClassName() << std::endl;
        find<SoftwareModel>()->run();
    }
};

#endif //SATURATIONPRESSUREMATERIALRELATION_H
