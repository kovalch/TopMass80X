#include <iostream>
#include <sstream> 
#include <boost/program_options.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/progress.hpp>
#include <vector>
#include <time.h>

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
    TString fLepton;
    TString fFileElectron;
    TString fFileMuon;
    TString fMethod;
    std::string fWeight;
    
    int fBins;
    double fLumi;
    double fSig;
    
    int target, run, luminosityBlock, event, combi, nVertex, leptonId;
    double hadTopMass, hadWRawMass, leptonPt, leptonC, hitFitProb, deltaThetaHadWHadB, deltaThetaHadQHadQBar, PUWeight, PUWeightUp, PUWeightDown, muWeight, bWeight, bWeight_bTagSFUp, bWeight_bTagSFDown, bWeight_misTagSFUp, bWeight_misTagSFDown, MCWeight, hadQBCSV, hadQBarBCSV, hadBBCSV, lepBBCSV;
    double jetsPt[4];
    double pdfWeights[44];

    TChain* fChain;
    
    TTree* fTree;
    TTree* fTreeTTmu;
    TTree* fTreeTTe;
    TTree* fTreeWmu;
    TTree* fTreeWe;
    TTree* fTreeSTmu;
    TTree* fTreeSTe;
    
    TFile* tempFile;
    TH2F* hEntries;
    TH2F* hMass;
    TH2F* hMassError;
    TH2F* hMassSigma;
    TH2F* hMassAlt;
    TH2F* hMassAltError;
    
    TH2F* hJES;
    TH2F* hJESError;
    
    TH2F* hMassCalibrated;
    TH2F* hMassErrorCalibrated;
    
    void CreateHistos();
    void CreateRandomSubset();
    void DrawEvents(TTree* tempTree, double nEventsPE);
    TTree* PrepareTree(TString file);

  public:
    Analysis(po::variables_map vm);
    Analysis(TString identifier, TString file, TString method, int bins, double lumi);
    ~Analysis() {}
    
    void Analyze(po::variables_map vm);
    
    TH2F* GetH2Mass();
    TH2F* GetH2MassError();
    TH2F* GetH2MassSigma();
    TH2F* GetH2MassAlt();
    TH2F* GetH2MassAltError();
    
    TH2F* GetH2JES();
    TH2F* GetH2JESError();
    
    TH2F* GetH2MassCalibrated();
    TH2F* GetH2MassErrorCalibrated();
    
    TString GetIdentifier();
};

