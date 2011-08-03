#include "MassAnalyzer.h"
#include "IdeogramCombLikelihood.h"

#include <iomanip>

#include "TH1D.h"
#include "TH2D.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TStyle.h"

class IdeogramAnalyzer : public MassAnalyzer {
  private:
    double QBTagProbability(double bDiscriminator);
  
  public:
    IdeogramAnalyzer(TString identifier, TTree* tree) : MassAnalyzer(identifier, tree) {};
    void Analyze(TString cuts, int i, int j);
    
    double GetMass();
};
