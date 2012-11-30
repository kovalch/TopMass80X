#include <iostream>
#include <sstream>

#include "TH1F.h"
#include "TString.h"
#include "TTree.h"

class Analysis {
 private:

  TString samplePath_;
  TString fIdentifier_;
  TString fMethod_;
  TString fBinning_;
  std::vector<float> vBinning_;

  //int fBins_;

  TString fChannel_;

  TTree* fTree_;

  std::map<TString, TH1F*> histograms_;

  void SetH1(TString histName, TH1F* hist);
  void CreateHisto(TString name);

 public:
  Analysis(std::vector<float> v);
  ~Analysis();
    
  void Analyze();

  TH1F* GetH1(TString histName);
  const std::map<TString, TH1F*> GetH1s() const;

  TString GetIdentifier();

};

