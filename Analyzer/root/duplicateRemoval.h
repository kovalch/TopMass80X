#ifndef DUPREMOVAL_C
#define DUPREMOVAL_C

#include "TTree.h"
#include "TChain.h"
#include "TH2F.h"

#include "RooFormulaVar.h"

#include <string>

class DuplicateRemover {
public:
  DuplicateRemover();
  ~DuplicateRemover() {}

private:

  void removeDuplicates();
};

#endif /* DUPREMOVAL_C */
