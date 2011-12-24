#include "MassAnalyzer.h"
#include "IdeogramCombLikelihood.h"
#include "Helper.h"

#include <iomanip>

#include "TH1D.h"
#include "TH2D.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TDirectory.h"

class IdeogramAnalyzer : public MassAnalyzer {
  private:
    double QBTagProbability(double bDiscriminator);
    double mWnVertex();
    double mTnVertex();
  
  public:
    IdeogramAnalyzer(TString identifier, TTree* tree) : MassAnalyzer(identifier, tree) {};
    void Analyze(TString cuts, int i, int j);
    
    double GetMass();
};
