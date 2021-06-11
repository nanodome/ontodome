#ifndef SURFACETENSIONPOLYNOMIALSOFTWAREMODEL_H
#define SURFACETENSIONPOLYNOMIALSOFTWAREMODEL_H

#include "../surfacetensionmaterialrelation.h"

class SurfaceTensionPolynomialSoftwareModel : public SoftwareModel {
private:
    /// Core implementation of the model.
    /// It stands fot the numerical representation of the model.
    /// \param T species temperature [K].
    double impl(double T)
    {
        return (s[0] - (s[1]*(T - s[2])));
    }

    std::vector<double> s; ///< Model coefficients array container.
    double* sval; ///< Pointer to the species's property data type.

public:
    SurfaceTensionPolynomialSoftwareModel() : SoftwareModel() {

//      SurfaceTension* dummyS = nullptr;
//      IUPAC* dummyI = nullptr;
//      Real* dummyR = nullptr;

//      this->SoftwareModel::createRelationTo<hasModel>(dummyS);
//      this->SoftwareModel::createRelationTo<hasInput>(dummyI);
//      this->SoftwareModel::createRelationTo<hasOutput>(dummyR);
    }

    std::string getClassName() const { return "SurfaceTensionPolynomialModel"; }

    /// Returns a brief description of the software model by means of its mathematical representation
    /// and the used parameters.
    std::string model_description() {
      return "a - b * ( T - c ) \nwhere a,b and c are the coefficients and T is the temperature.";
    }

    /// Select the model coefficients based on Species Symbol and quantity pointer.
    /// In order to save computational time on several calls of the run method, this
    /// is done on first call only.
    /// Then gathers all the necessary inputs through the relations graph and
    /// runs the actual implementation ("impl").
    void run() {
      if (s.size() == 0) {
        std::string _name = findNearest<SingleComponentComposition>()->name;
        s = get_coeffs(_name);
        sval = &findNearest<SurfaceTension>()->getRelatedObjects<Real>()[0]->data;
      }

      //Get the object's first related temperature
      double T = findNearest<Temperature>()->getRelatedObjects<Real>()[0]->data;

      // Compute the value and push it to the Species' Surface Tension object
      *sval = impl(T);
    }

    /// Return the coefficients needed by the model on given chemical element name.
    /// \param name species name.
    std::vector<double> get_coeffs(std::string _name) const {
      if (_name == "Silicon") {
        return {0.732, 0.000086, 1685};
      }
      else if (_name == "Argon") {
        return {0.,0.,0.};
      }
      else if (_name == "Helium") {
        return {0.,0.,0.};
      }
      else { abort(); }
    }

};

#endif // SURFACETENSIONPOLYNOMIALSOFTWAREMODEL_H
