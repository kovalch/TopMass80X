#include <iostream>
#include <sstream>

#include "TH2F.h"
#include "TString.h"
#include "TTree.h"

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
  Analysis(std::vector<float> v);
  ~Analysis();
    
  void Analyze();

  TH1F* GetH1(TString histName);
  const std::map<TString, TH2F*> GetH2s() const;

  TString GetIdentifier();

};

