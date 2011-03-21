#include "MassAnalyzer.h"
#include "IdeogramCombLikelihood.h"

#include "TF1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TMath.h"

class IdeogramAnalyzer : public MassAnalyzer {
  private:
  
  public:
    IdeogramAnalyzer(TString identifier, TChain* chain) : MassAnalyzer(identifier, chain) {};
    void Analyze(TString cuts, int i, int j);
    
    double GetMass();
};
