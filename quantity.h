#ifndef QUANTITY_H
#define QUANTITY_H

#include <valarray>

#include "baseclass.h"

class Quantity : public BaseClass
{
    std::string unit;
public:
    Quantity(std::string _unit);

    std::string getClassName() { return "Quantity"; }

    std::string getUnit() const;
};

class ScalarQuantity : public Quantity
{
    double value;
public:
    ScalarQuantity(double _value, std::string _unit);

    std::string getClassName() { return "ScalarQuantity"; }

    double getValue() const;
};

class VectorQuantity : public Quantity
{
    std::valarray<double> value;
public:
    VectorQuantity(const std::valarray<double>& _value, std::string _unit);

    std::string getClassName() { return "VectorQuantity"; }

    std::valarray<double> getValue() const;
};

#endif // QUANTITY_H
