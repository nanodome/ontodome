#ifndef SATURATIONPRESSUREMODEL_H
#define SATURATIONPRESSUREMODEL_H

#include "../../base/thing.h"

class SaturationPressureModel : public SoftwareModel {
protected:
    double impl(std::vector<double> s, double T)
    {
        return pow(10.0,(s[0]-(s[1]/T)));
    }

    std::vector<double> s;

public:
    SaturationPressureModel(std::vector<double> _s) : SoftwareModel() { s = _s; }

    std::string getClassName() const { return "SaturationPressureModel"; }

    void run() {

      auto spec = this->getRelations<hasModel>()[0]->getDomain()->getRelations<hasProperty>()[0]->getDomain();

      double T = spec->getRelatedScalarObjects<Temperature>()[0];

      spec->getRelatedObjects<SaturationPressure>()[0]->getRelatedObjects<Scalar>()[0]->data = impl(s,T);
    }
};

#endif // SATURATIONPRESSUREMODEL_H
