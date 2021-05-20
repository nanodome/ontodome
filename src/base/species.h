#include "thing.h"

class HeliumSymbol : public ChemicalElement {
  static std::string name;

public:
  static std::string get_symbol() { return name; }
};

std::string HeliumSymbol::name = "Helium";

class Helium : public Atom {
public:
    Helium () : Atom () {
      HeliumSymbol* he = 0;
      this->createRelationTo<hasProperty>(he);
    }
    std::string getClassName() const { return "Helium"; }
};

class ArgonSymbol : public ChemicalElement {
  static std::string name;

public:
  static std::string get_symbol() { return name; }
};

std::string ArgonSymbol::name = "Argon";

class Argon : public Atom {
public:
    Argon () : Atom () {
      ArgonSymbol* ar = 0;
      this->createRelationTo<hasProperty>(ar);
    }
    std::string getClassName() const { return "Argon"; }
};

class SiliconSymbol : public ChemicalElement {
  static std::string name;

public:
  static std::string get_symbol() { return name; }
};

std::string SiliconSymbol::name = "Silicon";

class Silicon : public Atom {
public:
    Silicon () : Atom () {
      SiliconSymbol* si = 0;
      this->createRelationTo<hasProperty>(si);
    }
    std::string getClassName() const { return "Silicon"; }
};
