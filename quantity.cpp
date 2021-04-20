#include "quantity.h"

Quantity::Quantity(std::string _unit) : unit(_unit) {}

std::string Quantity::getUnit() const { return unit; }

ScalarQuantity::ScalarQuantity(double _value, std::string _unit) :
    Quantity(_unit), value(_value) {}

double ScalarQuantity::getValue() const { return value; }


VectorQuantity::VectorQuantity(const std::valarray<double>& _value, std::string _unit) :
    Quantity(_unit), value(_value) {}

std::valarray<double> VectorQuantity::getValue() const { return value; }
