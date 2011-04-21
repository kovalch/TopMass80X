#include "MassAnalyzer.h"
#include "IdeogramCombLikelihood.h"

#include "TH1D.h"
#include "TCanvas.h"
#include "TMath.h"

class IdeogramAnalyzer : public MassAnalyzer {
  private:
    double QBTagProbability(double bDiscriminator);
  
  public:
    IdeogramAnalyzer(TString identifier, TTree* tree) : MassAnalyzer(identifier, tree) {};
    void Analyze(TString cuts, int i, int j);
    
    double GetMass();
};
