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
#include "IdeogramAnalyzer.h"

class Analysis {
  private:
    TString fIdentifier;
    TString fFile;
    TString fMethod;
    MassAnalyzer* fAnalyzer;
    int fBins;
    
    bool fAnalyzed;
    
    TChain* fChain;
    TH2F* hEntries;
    TH2F* hMass;
    TH2F* hMassError;
    TH2F* hMassSigma;
    
    TH2F* hMassCalibrated;
    TH2F* hMassErrorCalibrated;
    
    void CreateHistos();
    TTree* CreateRandomSubset();

  public:
    Analysis(TString identifier, TString file, TString method, int bins) :
      fIdentifier(identifier), fFile(file), fMethod(method), fBins(bins) {}
    ~Analysis() { delete &fAnalyzed; }
    
    void Analyze(bool reanalyze = false);
    
    TH2F* GetH2Mass();
    TH2F* GetH2MassError();
    TH2F* GetH2MassSigma();
    
    TH2F* GetH2MassCalibrated();
    TH2F* GetH2MassErrorCalibrated();
    
    TString GetIdentifier();
};

