#include <iostream>
#include <sstream> 

#include "TCanvas.h"
#include "TChain.h"
#include "TH2F.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"

#include "GenMatchAnalyzer.h"
#include "MVAAnalyzer.h"

class Analysis {
  private:
    TString fIdentifier;
    TString fFile;
    TString fMethod;
    MassAnalyzer* fAnalyzer;
    int fBins;
    
    TChain* fChain;
    TH2F* hEntries;
    TH2F* hMass;
    TH2F* hMassError;
    TH2F* hMassSigma;
    
    void CreateHistos();

  public:
    Analysis(TString identifier, TString file, TString method, int bins) :
      fIdentifier(identifier), fFile(file), fMethod(method), fBins(bins) {}
    void Analyze();
    
    TH2F* GetH2Mass();
    TH2F* GetH2MassError();
    TH2F* GetH2MassSigma();
    
    TString GetIdentifier();
};

