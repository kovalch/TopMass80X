#include <vector>
#include <cmath>
#include <fstream>

#include "Analysis.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH3F.h"
#include "TVector.h"

bool fexists(const char *filename)
{
  ifstream ifile(filename);
  return ifile;
}

struct massPoint {
  double genMass;
  TString identifier;
  TString fileName;
  Analysis* analysis;
  
  TH3F* h3Mass;
  TH3F* h3MassPull;
  
  massPoint(double pGenMass, TString pIdentifier) :
      genMass(pGenMass), identifier(pIdentifier) {
    if (fexists("/scratch/hh/lustre/cms/user/mseidel/root/analyzeTop_1725.root")) {
      fileName = "/scratch/hh/lustre/cms/user/mseidel/root/analyzeTop_";
    }
    else fileName = "root/analyzeTop_";
    fileName += identifier;
    fileName += ".root";
  };
};

class TopMass {
  private:
    TString fMethod;
    int fBins;
    double fLumi;
    
    std::vector<massPoint> massPoints;
    std::vector<massPoint>::iterator iMassPoint;
    
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
