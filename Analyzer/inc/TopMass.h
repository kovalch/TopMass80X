#include <vector>
#include <cmath>

#include "Analysis.h"
#include "TGraphErrors.h"

class TopMass {
  private:
    TString fMethod;
    int fBins;
    double fLumi;
    
    std::vector<Analysis*> calibrationAnalyses;
    std::vector<Analysis*>::const_iterator iAnalysis;
    
    Analysis* a1665;
    Analysis* a1725;
    Analysis* a1785;
    
    Analysis* a1665_jes_up;
    Analysis* a1725_jes_up;
    Analysis* a1785_jes_up;
    
    Analysis* a1665_jes_down;
    Analysis* a1725_jes_down;
    Analysis* a1785_jes_down;
    
    Analysis* aSim;
    
    double fCalibFitParameter[6][6][2];
    double fCalibFitParError[6][6][2];
  
  public:
    TopMass(TString method, int bins, double lumi);
    
    void WriteEnsembleTestTree();
    void EvalEnsembleTest();
    void Calibrate();
    TH2F* Measure(Analysis* a);
    void Systematics();
};
