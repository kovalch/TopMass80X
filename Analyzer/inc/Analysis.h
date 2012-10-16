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
    TString fIdentifier, fLepton, fFileElectron, fFileMuon, fMethod, fBinning;
    std::string fWeight;
    
    int fBins;
    double fLumi;
    double fSig;
    std::vector<float> vBinning;
    
    int target, run, luminosityBlock, event, combi, nVertex, leptonId;
    double hadTopMass, hadWRawMass, leptonPt, leptonC, hitFitProb, deltaThetaHadWHadB, deltaThetaHadQHadQBar, PUWeight, PUWeightUp, PUWeightDown, muWeight, bWeight, bWeight_bTagSFUp, bWeight_bTagSFDown, bWeight_misTagSFUp, bWeight_misTagSFDown, MCWeight, mcWeight, hadQBCSV, hadQBarBCSV, hadBBCSV, lepBBCSV;
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
    TH1F* hEntries;
    TH1F* hMass;
    TH1F* hMassError;
    TH1F* hMassSigma;
    TH1F* hMassAlt;
    TH1F* hMassAltError;
    
    TH1F* hJES;
    TH1F* hJESError;
    
    TH1F* hMassCalibrated;
    TH1F* hMassErrorCalibrated;
    
    void CreateHistos();
    void CreateRandomSubset();
    void DrawEvents(TTree* tempTree, double nEventsPE);
    TTree* PrepareTree(TString file);

  public:
    Analysis(po::variables_map vm, std::vector<float> v);
    Analysis(TString identifier, TString file, TString method, int bins, double lumi);
    ~Analysis() {}
    
    void Analyze(po::variables_map vm);
    
    TH1F* GetH1Mass();
    TH1F* GetH1MassError();
    TH1F* GetH1MassSigma();
    TH1F* GetH1MassAlt();
    TH1F* GetH1MassAltError();
    
    TH1F* GetH1JES();
    TH1F* GetH1JESError();
    
    TH1F* GetH1MassCalibrated();
    TH1F* GetH1MassErrorCalibrated();
    
    TString GetIdentifier();
};

