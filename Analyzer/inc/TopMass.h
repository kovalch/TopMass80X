#ifndef TOPMASS_H
#define TOPMASS_H

#include <vector>

#include "TString.h"

class TopMass {
  private:
    //bool fexists(const char *filename);
    void WriteEnsembleTest(const std::vector<float>& vBinning);
      
    const TString fBinning_, fTask_;

  public:
    TopMass();
};

#endif
