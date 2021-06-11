#ifndef SATURATIONPRESSUREPOLYNOMIALSOFTWAREMODEL_H
#define SATURATIONPRESSUREPOLYNOMIALSOFTWAREMODEL_H

#include "../saturationpressurematerialrelation.h"

class SaturationPressurePolynomialSoftwareModel : public SoftwareModel {
private:
    double impl(double T)
    {
        return 1.01e5 * pow(10.0,(s[0]-(s[1]/T)));
    }

    std::vector<double> s;
    double* sval;

public:
    SaturationPressurePolynomialSoftwareModel() : SoftwareModel() {

//      SaturationPressure* dummyS = nullptr;
//      IUPAC* dummyI = nullptr;
//      Real* dummyR = nullptr;

//      this->SoftwareModel::createRelationTo<hasModel>(dummyS);
//      this->SoftwareModel::createRelationTo<hasInput>(dummyI);
//      this->SoftwareModel::createRelationTo<hasOutput>(dummyR);
    }

    std::string getClassName() const { return "SaturationPressurePolynomialModel"; }

    std::string model_description() {
      return "10^(a - b/T) \nwhere a and b are the coefficients and T is the temperature.";
    }

    void run() {
      // Select the model coefficients based on Species Symbol and quantity pointer
      // In order to save computational time on several calls of the run method, this
      // has to be done the first time only
      if (s.size() == 0) {
        std::string _name = findNearest<SingleComponentComposition>()->name;
        s = get_coeffs(_name);
        sval = &findNearest<SaturationPressure>()->getRelatedObjects<Real>()[0]->data;
      }

      //Get the object's first related temperature
      double T = findNearest<Temperature>()->getRelatedObjects<Real>()[0]->data;

      // Compute the value and push it to the Species' Saturation Pressure object
      *sval = impl(T);
    }

    std::vector<double> get_coeffs(std::string _name) const {
      if (_name == "Silicon") {
        return {7.5341, 23399.};
      }
      else if (_name == "Argon") {
        return {0.,0.};
      }
      else if (_name == "Helium") {
        return {0.,0.};
      }
      else { abort(); }
    }
};

#endif // SATURATIONPRESSUREPOLYNOMIALSOFTWAREMODEL_H
