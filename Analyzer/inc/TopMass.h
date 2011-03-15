#include <vector>
#include <cmath>

#include "Analysis.h"
#include "TGraphErrors.h"

class TopMass {
  private:
    int fBins;
    
    std::vector<Analysis*> analyses;
    std::vector<Analysis*>::const_iterator iAnalysis;
    
    double fCalibFitParameter[8][8][2];
    double fCalibFitParError[8][8][2];
  
  public:
    TopMass(int bins) : fBins(bins) {};
    void Calibrate(TString method);
    void Measure(TString method);
};
