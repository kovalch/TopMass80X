#ifndef MASSANALYZER_H
#define MASSANALYZER_H

#include "TTree.h"
#include "TF1.h"

class MassAnalyzer {
  protected:
    TString fIdentifier;
    TTree* fTree;
    double fEntries;
    double fMass;
    double fMassError;
    double fMassSigma;
    double fMassAlt;
    double fMassAltError;
    double fJES;
    double fJESError;
  
  public:
    MassAnalyzer(TString identifier, TTree* tree) : fIdentifier(identifier), fTree(tree) {};
    ~MassAnalyzer() { delete fTree; }
    
    
    virtual void Analyze(TString cuts, int i, int j) = 0;
    
    double GetMass();
    double GetMassError();
    double GetMassSigma();
    double GetMassAlt();
    double GetMassAltError();
    double GetJES();
    double GetJESError();
};

#endif
