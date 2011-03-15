#ifndef MASSANALYZER_H
#define MASSANALYZER_H

#include "TChain.h"
#include "TF1.h"

class MassAnalyzer {
  protected:
    TString fIdentifier;
    TChain* fChain;
    double fEntries;
    double fMass;
    double fMassError;
    double fMassSigma;
  
  public:
    MassAnalyzer(TString identifier, TChain* chain) : fIdentifier(identifier), fChain(chain) {};
    
    virtual void Analyze(TString cuts, int i, int j) = 0;
    
    double GetMass();
    double GetMassError();
    double GetMassSigma();
};

#endif
