#ifndef INFREMOVAL_C
#define INFREMOVAL_C

#include "RooFormulaVar.h"

#include <string>

class InfRemover {
public:
  InfRemover();
  ~InfRemover() {}

private:

  void removeInfinities();
};

#endif /* INFREMOVAL_C */
