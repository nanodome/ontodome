#include "quantity.h"

Quantity::Quantity(std::string _unit, std::string _name) : unit(_unit), name(_name) {}

std::string Quantity::getUnit() const { return unit; }

std::string Quantity::getName() const { return name; }

ScalarQuantity::ScalarQuantity(double _value, std::string _unit, std::string _name) :
    Quantity(_unit, _name), value(_value) {}

double ScalarQuantity::getValue() { return value; }

VectorQuantity::VectorQuantity(const std::valarray<double> _value, std::string _unit, std::string _name) :
    Quantity(_unit, _name), value(_value) {}

std::valarray<double> VectorQuantity::getValue() { return value; }
