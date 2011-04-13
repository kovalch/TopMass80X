#include <iostream>
#include <sstream> 

#include "TCanvas.h"
#include "TChain.h"
#include "TH2F.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TFile.h"

#include "Helper.h"
#include "GenMatchAnalyzer.h"
#include "MVAAnalyzer.h"
#include "IdeogramAnalyzer.h"

class Analysis {
  private:
    TString fIdentifier;
    TString fFile;
    TString fMethod;

    int fBins;
    double fLumi;
    
    bool fAnalyzed;
    
    TChain* fChain;
    TTree* fTree;
    TH2F* hEntries;
    TH2F* hMass;
    TH2F* hMassError;
    TH2F* hMassSigma;
    
    TH2F* hMassCalibrated;
    TH2F* hMassErrorCalibrated;
    
    void CreateHistos();
    void CreateRandomSubset();

  public:
    Analysis(TString identifier, TString file, TString method, int bins, double lumi);
    ~Analysis() { delete &fAnalyzed; }
    
    void Analyze(bool reanalyze = false);
    
    TH2F* GetH2Mass();
    TH2F* GetH2MassError();
    TH2F* GetH2MassSigma();
    
    TH2F* GetH2MassCalibrated();
    TH2F* GetH2MassErrorCalibrated();
    
    TString GetIdentifier();
};

