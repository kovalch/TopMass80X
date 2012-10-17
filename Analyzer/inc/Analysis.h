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

#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "Helper.h"
#include "GenMatchAnalyzer.h"
#include "MVAAnalyzer.h"
#include "IdeogramAnalyzer.h"
#include "RooFitTemplateAnalyzer.h"

namespace po = boost::program_options;

class Analysis {
 private:
  TString samplePath;
  TString fIdentifier;
  TString fFile;
  TString fMethod;

  int fBins;
  double fLumi;
    
  TChain* fChain;
  TTree* fTree;
  TTree* tTree;
  TTree* tTreeBkg;
  TFile* tempFile;
  TH2F* hEntries;

  TH2F* hMass;
  TH2F* hMassError;
  //TH2F* hMassSigma;
  TH2F* hJES;
  TH2F* hJESError;
    
  TH2F* hMassConstJES;
  TH2F* hMassConstJESError;
  //TH2F* hMassConstJESSigma;

  TH2F* hFSig;
  TH2F* hFSigError;
  TH2F* hMassfSig;
  TH2F* hMassfSigError;
  //TH2F* hMassfSigSigma;
  TH2F* hJESfSig;
  TH2F* hJESfSigError;

  //TH1F * bTagEffScaleFactor;
  
  TH2F * bTagEff;
  TH2F * cTagEff;
  TH2F * lTagEff;

  //TH2F* hMassCalibrated;
  //TH2F* hMassErrorCalibrated;
    
  void CreateHistos();
  void CreateRandomSubset();

 public:
  Analysis(po::variables_map vm);
  Analysis(TString identifier, TString file, TString method, int bins, double lumi);
  ~Analysis();
    
  void Analyze(po::variables_map vm);

  //enum enumForPUWeights {kSummer11, kSummer11Plus05, kSummer11Minus05, kFall11, kFall11Plus05, kFall11Minus05, kFall10};
  void AddWeights  (TTree* tempTree, bool isData=false); //, enumForPUWeights whichSample, int whichPDF = 0);
  void AdaptBkgTree(TTree* tempTreeBkg);
    
  TH2F* GetH2Mass();
  TH2F* GetH2MassError();
  //TH2F* GetH2MassSigma();
  TH2F* GetH2JES();
  TH2F* GetH2JESError();
    
  TH2F* GetH2MassConstJES();
  TH2F* GetH2MassConstJESError();
  //TH2F* GetH2MassConstJESSigma();
    
  TH2F* GetH2FSig();
  TH2F* GetH2FSigError();
  TH2F* GetH2MassfSig();
  TH2F* GetH2MassfSigError();
  //TH2F* GetH2MassfSigSigma();
  TH2F* GetH2JESfSig();
  TH2F* GetH2JESfSigError();
    
  //TH2F* GetH2MassCalibrated();
  //TH2F* GetH2MassErrorCalibrated();
    
  TString GetIdentifier();

  // return the PU weights for the different samples
  enum enumForPUWeights {kSummer11, kSummer11Plus05, kSummer11Minus05, kFall11, kFall11Plus05, kFall11Minus05, kFall10};
  double calcPUWeight_(enum enumForPUWeights sample, short nPU);
  double calcPDFWeight_(int whichPDFUncertainty, bool upVariation, double x1, int id1, float Q, double x2, int id2);
  double eventBTagProbability_(std::vector<double> &oneMinusBEffies, std::vector<double> &oneMinusBMistags, bool verbose = false);
  double calcBTagWeight_(int Njet=0, short * pdgId=0, TClonesArray * jets=0);
};

