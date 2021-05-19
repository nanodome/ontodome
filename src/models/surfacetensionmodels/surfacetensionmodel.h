#ifndef SURFACETENSIONMODEL_H
#define SURFACETENSIONMODEL_H

#include "../../base/thing.h"

class SurfaceTensionModel : public SoftwareModel {
protected:
    double model(std::vector<double> s, double T)
    {
        return (s[0] - (s[1]*(T - s[2])));
    }

    std::vector<double> s;

public:
    // Dummy constructor for Knowledge Generators navigation
    SurfaceTensionModel() : SoftwareModel() {
      this->createRelationTo<isModelFor>(new GasMixture);
      this->createRelationTo<isModelFor>(new SurfaceTension(new Scalar(0), new Unit("dummy")));
      this->createRelationTo<hasInput>(new Vector({0}));
      this->createRelationTo<hasOutput>(new Scalar(0));

      this->createRelationTo<isSoftwareModelFor,ContinuumModel>(new ContinuumModel(
            "a - b * ( T - c ) \nwhere a,b and c are the coefficients and T is the temperature.")
        );
    };

    // This model requires an initialization which is a vector containing the model's coefficients
    SurfaceTensionModel(std::vector<double> _s) : SoftwareModel() {

      // Store the model coefficients
      s = _s;

      this->createRelationTo<isModelFor>(new GasMixture); // actual relation
      this->createRelationTo<isModelFor>(new SurfaceTension(new Scalar(0), new Unit("dummy")));
      this->createRelationTo<hasInput>(new Vector({0}));
      this->createRelationTo<hasOutput>(new Scalar(0));

      this->createRelationTo<isSoftwareModelFor,ContinuumModel>(new ContinuumModel(
            "a - b * ( T - c ) \nwhere a,b and c are the coefficients and T is the temperature.")
        );
    }

    std::string getClassName() const { return "SurfaceTensionModel"; }

    void run(double T) {

      // Compute the value and push it to the Species' Saturation pressure property
      this->getLastRelation<hasModel>()->getDomain()->getRelatedObjects<SurfaceTension>()[0]->getRelatedObjects<Scalar>()[0]->data = model(s,T);
    }
};

#endif // SURFACETENSIONMODEL_H
