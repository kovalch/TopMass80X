#include <iostream>
#include <sstream> 

#include "TCanvas.h"
#include "TChain.h"
#include "TH2F.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"

#include "DummyAnalyzer.h"
#include "GenMatchAnalyzer.h"

class Analysis {
  private:
    int fBins;
    TChain* fChain;
    TH2F* hEntries;
    TH2F* hMassMean;
    TH2F* hMassMeanError;
    TH2F* hMassSigma;
    
    void CreateHistos();

  public:
    Analysis() { fBins = 8; }
    void Print() const;
    void SetBins(int x) { fBins = x; }
    
    void Analyze();
};

