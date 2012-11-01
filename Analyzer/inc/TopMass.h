#ifndef TOPMASS_H
#define TOPMASS_H

#include <vector>

#include "TString.h"

#include "Analysis.h"

class TopMass {
  private:
    TString fMethod_;
    TString fBinning_;
    double fLumi_;

    bool fexists(const char *filename);
      
  public:
    TopMass();
    
    void WriteEnsembleTest(std::vector<float> vBinning);
};

#endif
