#ifndef MASSANALYZER_H
#define MASSANALYZER_H

#include "TChain.h"
#include "TF1.h"

class MassAnalyzer {
  protected:
    TChain* fChain;
    double fEntries;
    double fMass;
    double fMassError;
    double fMassSigma;
  
  public:
    MassAnalyzer(TChain* chain) { fChain = chain; };
    
    virtual void Analyze(TString cuts, int i, int j) = 0;
    
    double GetMass();
    double GetMassError();
    double GetMassSigma();
};

#endif
