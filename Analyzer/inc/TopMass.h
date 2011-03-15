#include <vector>
#include <cmath>

#include "Analysis.h"
#include "TGraphErrors.h"

class TopMass {
  private:
    TString fMethod;
    int fBins;
    
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
    
    double fCalibFitParameter[8][8][2];
    double fCalibFitParError[8][8][2];
  
  public:
    TopMass(TString method, int bins);
    
    void Calibrate();
    TH2F* Measure(Analysis* a);
    void Systematics();
};
