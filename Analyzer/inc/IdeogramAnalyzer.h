#include "MassAnalyzer.h"
#include "TF1.h"
#include "TCanvas.h"

class IdeogramAnalyzer : public MassAnalyzer {
  private:
    TF1* gaus;
  
  public:
    IdeogramAnalyzer(TString identifier, TChain* chain) : MassAnalyzer(identifier, chain) {};
    void Analyze(TString cuts, int i, int j);
    
    double GetMass();
};
