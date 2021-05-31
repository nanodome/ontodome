#ifndef SURFACETENSIONPOLYNOMIALMODEL_H
#define SURFACETENSIONPOLYNOMIALMODEL_H

#include "../surfacetensionmaterialrelation.h"

class SurfaceTensionPolynomialModel : public SoftwareModel, public SurfaceTensionMaterialRelation {
private:
    double impl(double T)
    {
        return (s[0] - (s[1]*(T - s[2])));
    }

    std::vector<double> s;

public:
    SurfaceTensionPolynomialModel() : SoftwareModel() {

//      SurfaceTension* dummyS = nullptr;
//      IUPAC* dummyI = nullptr;
//      Real* dummyR = nullptr;

//      this->SoftwareModel::createRelationTo<hasModel>(dummyS);
//      this->SoftwareModel::createRelationTo<hasInput>(dummyI);
//      this->SoftwareModel::createRelationTo<hasOutput>(dummyR);
    }

    std::string getClassName() const { return "SurfaceTensionPolynomialModel"; }

    std::string model_description() {
      return "a - b * ( T - c ) \nwhere a,b and c are the coefficients and T is the temperature.";
    }

    void run() {
      // Select the model coefficients based on IUPAC name
//      std::string _name = "Si";
      std::string _name = find<SingleComponentComposition>()->name;
      s = get_coeffs(_name);

      //Get the object's first related temperature
      double T = find<Temperature>()->getRelatedObjects<Real>()[0]->data;

      // Compute the value and push it to the Species' Surface Tension object
      find<SurfaceTension>()->getRelatedObjects<Real>()[0]->data = impl(T);
    }

    std::vector<double> get_coeffs(std::string name) const {
      if (name == "Silicon") {
        return {0.732, 0.000086, 1685};
      }
      else if (name == "Argon") {
        return {0.,0.,0.};
      }
      else if (name == "Helium") {
        return {0.,0.,0.};
      }
      else { abort(); }
    }
};

#endif // SURFACETENSIONPOLYNOMIALMODEL_H
