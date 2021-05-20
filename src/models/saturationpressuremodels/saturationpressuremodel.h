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
      this->createRelationTo<hasModel>(new SaturationPressure(new Scalar(0), new Unit("dummy")));
      this->createRelationTo<hasInput>(new IUPAC(""));
      this->createRelationTo<hasInput>(new Scalar(0));
      this->createRelationTo<hasOutput>(new Scalar(0));

      this->createRelationTo<hasSoftwareModel,ContinuumModel>(new ContinuumModel(
            "10^(a - b/T) \nwhere a and b are the coefficients and T is the temperature.")
            );
    };

    // This model requires an initialization which is a vector containing the model's coefficients
//    SaturationPressureModel(std::vector<double> _s) : SoftwareModel() {
      SaturationPressureModel(ChemicalSpecies* chem) : SoftwareModel() {

      std::string name = chem->getRelatedObjects<IUPAC>()[0];
      // Store the model coefficients
//      std::vector<double> s = get_coeffs(_s);
//      s = _s;
      s = {0.1,0.2};

      this->createRelationTo<hasModel>(new SaturationPressure(new Scalar(0), new Unit("dummy")));
      this->createRelationTo<hasInput>(new IUPAC(""));
      this->createRelationTo<hasInput>(new Scalar(0));
      this->createRelationTo<hasOutput>(new Scalar(0));

      this->createRelationTo<hasSoftwareModel,ContinuumModel>(new ContinuumModel(
            "10^(a - b/T) \nwhere a and b are the coefficients and T is the temperature.")
            );
    }

    std::string getClassName() const { return "SaturationPressureModel"; }

    void run(double T) {

      // Compute the value and push it to the Species' Saturation pressure property
      this->getLastRelation<hasModel>()->getDomain()->getRelatedObjects<SaturationPressure>()[0]->getRelatedObjects<Scalar>()[0]->data = model(s,T);
    }
};

#endif // SATURATIONPRESSUREMODEL_H
