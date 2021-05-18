#ifndef SATURATIONPRESSUREMODEL_H
#define SATURATIONPRESSUREMODEL_H

#include "../../base/thing.h"

class SaturationPressureModel : public SoftwareModel {
protected:
    double model(std::vector<double> s, double T)
    {
        return pow(10.0,(s[0]-(s[1]/T)));
    }

    std::vector<double> s;

public:
    // Dummy constructor for Knowledge Generators navigation
    SaturationPressureModel() : SoftwareModel() {
      this->createRelationTo<isModelFor>(new GasMixture);
      this->createRelationTo<isModelFor>(new SaturationPressure(new Scalar(0), new Unit("dummy")));
    };

    // This model requires an initialization which is a vector containing the model's coefficients
    SaturationPressureModel(std::vector<double> _s) : SoftwareModel() {

      // Store the model coefficients
      s = _s;

      this->createRelationTo<isModelFor>(new GasMixture);
      this->createRelationTo<isModelFor>(new SaturationPressure(new Scalar(0), new Unit("dummy")));
    }

    std::string getClassName() const { return "SaturationPressureModel"; }

    void run(double T) {

      // Compute the value and push it to the Species' Saturation pressure property
      this->getLastRelation<hasModel>()->getDomain()->getRelatedObjects<SaturationPressure>()[0]->getRelatedObjects<Scalar>()[0]->data = model(s,T);
    }
};

#endif // SATURATIONPRESSUREMODEL_H
