#include "MassAnalyzer.h"
#include "TF1.h"
#include "TCanvas.h"

class GenMatchAnalyzer : public MassAnalyzer {
  private:
    TF1* gaus;
  
  public:
    GenMatchAnalyzer(TString identifier, TTree* tree) : MassAnalyzer(identifier, tree) {};
    void Analyze(TString cuts, int i, int j);
    
    double GetMass();
};
