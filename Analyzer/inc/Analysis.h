#include <iostream>
#include <sstream> 
#include <boost/program_options.hpp>

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

namespace po = boost::program_options;

class Analysis {
  private:
    TString fIdentifier;
    TString fFile;
    TString fMethod;

    int fBins;
    double fLumi;
    
    TChain* fChain;
    TTree* fTree;
    TTree* tTree;
    TFile* tempFile;
    TH2F* hEntries;
    TH2F* hMass;
    TH2F* hMassError;
    TH2F* hMassSigma;
    
    TH2F* hJES;
    TH2F* hJESError;
    
    TH2F* hMassCalibrated;
    TH2F* hMassErrorCalibrated;
    
    void CreateHistos();
    void CreateRandomSubset();

  public:
    Analysis(po::variables_map vm);
    Analysis(TString identifier, TString file, TString method, int bins, double lumi);
    ~Analysis() {}
    
    void Analyze(po::variables_map vm);
    
    TH2F* GetH2Mass();
    TH2F* GetH2MassError();
    TH2F* GetH2MassSigma();
    
    TH2F* GetH2JES();
    TH2F* GetH2JESError();
    
    TH2F* GetH2MassCalibrated();
    TH2F* GetH2MassErrorCalibrated();
    
    TString GetIdentifier();
};

