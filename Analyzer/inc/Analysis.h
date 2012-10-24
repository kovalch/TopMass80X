#include <iostream>
#include <sstream>
#include "boost/program_options.hpp"

#include "TCanvas.h"
#include "TChain.h"
#include "TH2F.h"
#include "TF1.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TFile.h"

#include "TClonesArray.h"
#include "TLorentzVector.h"

namespace po = boost::program_options;

class Analysis {
 private:
  TString _samplePath;
  TString _fIdentifier;
  TString _fMethod;

  int _fBins;

  TString _fChannel;

  TTree* _fTree;
  TH2F* _hEntries;

  TH2F* _hMass;
  TH2F* _hMassError;
  TH2F* _hJES;
  TH2F* _hJESError;
    
  TH2F* _hMassConstJES;
  TH2F* _hMassConstJESError;

  TH2F* _hFSig;
  TH2F* _hFSigError;
  TH2F* _hMassfSig;
  TH2F* _hMassfSigError;
  TH2F* _hJESfSig;
  TH2F* _hJESfSigError;

  void CreateHistos();

 public:
  Analysis(po::variables_map vm);
  Analysis(TString identifier, TString file, TString method, int bins, double lumi);
  ~Analysis();
    
  void Analyze(po::variables_map vm);

  TH2F* GetH2Mass();
  TH2F* GetH2MassError();
  TH2F* GetH2JES();
  TH2F* GetH2JESError();
    
  TH2F* GetH2MassConstJES();
  TH2F* GetH2MassConstJESError();
    
  TH2F* GetH2FSig();
  TH2F* GetH2FSigError();
  TH2F* GetH2MassfSig();
  TH2F* GetH2MassfSigError();
  TH2F* GetH2JESfSig();
  TH2F* GetH2JESfSigError();
    
  TString GetIdentifier();

};

