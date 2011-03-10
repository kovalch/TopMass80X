#include "MassAnalyzer.h"
#include "TF1.h"
#include "TCanvas.h"

class GenMatchAnalyzer : public MassAnalyzer {
  private:
    TF1* gaus;
  
  public:
    GenMatchAnalyzer(TChain* chain) : MassAnalyzer(chain) {};
    void Analyze(TString cuts);
};
