#ifndef SURFACETENSIONMODEL_H
#define SURFACETENSIONMODEL_H

#include "../../base/thing.h"

class SurfaceTensionModel : public SoftwareModel {
protected:
    double impl(std::vector<double> s, double T)
    {
        return (s[0] - (s[1]*(T - s[2])));
    }

    std::vector<double> s;

public:
    SurfaceTensionModel(std::vector<double> _s) : SoftwareModel() { s = _s; }

    std::string getClassName() const { return "SurfaceTensionModel"; }

    void run() {

        auto spec = this->getRelations<hasModel>()[0]->getDomain()->getRelations<hasProperty>()[0]->getDomain();

        double T = spec->getRelatedScalarObjects<Temperature>()[0];

        spec->getRelatedObjects<SurfaceTension>()[0]->getRelatedObjects<Scalar>()[0]->data = impl(s,T);
    }
};

#endif // SURFACETENSIONMODEL_H
