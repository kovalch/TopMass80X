#include <iostream>
#include <sstream>
#include "boost/program_options.hpp"

#include "TH2F.h"
#include "TString.h"
#include "TTree.h"

namespace po = boost::program_options;

class Analysis {
 private:
  TString samplePath_;
  TString fIdentifier_;
  TString fMethod_;

  int fBins_;

  TString fChannel_;

  TTree* fTree_;

  std::map<TString, TH2F*> histograms_;

  void SetH2(TString histName, TH2F* hist);
  void CreateHisto(TString name);

 public:
  Analysis(po::variables_map vm);
  Analysis(TString identifier, TString file, TString method, int bins, double lumi);
  ~Analysis();
    
  void Analyze(po::variables_map vm);

  TH2F* GetH2(TString histName);
  const std::map<TString, TH2F*> GetH2s() const;

  TString GetIdentifier();

};

