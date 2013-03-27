#ifndef TOPMASS_H
#define TOPMASS_H

#include <vector>

#include "TString.h"

#include "Analysis.h"

class TopMass {
  private:
    bool fexists(const char *filename);
    void WriteEnsembleTest(std::vector<float> vBinning);
      
    const TString fBinning_, fTask_;

  public:
    TopMass();
};

#endif
