#include <vector>
#include <cmath>
#include <fstream>

#include "Analysis.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH3F.h"
#include "TVector.h"

#include "tinyxml.h"

#include "Helper.h"

bool fexists(const char *filename)
{
  ifstream ifile(filename);
  return ifile;
}

struct massPoint {
  double genMass;
  double genLumi;
  TString identifier;
  TString fileName;
  Analysis* analysis;
  
  TH2F* h2Mass;
  TH2F* h2MassError;
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
    
    Analysis* aSim;
    
    double fCalibFitParameter[6][6][2];
    double fCalibFitParError[6][6][2];
    
    void LoadXML();
  
  public:
    TopMass(TString method, int bins, double lumi);
    
    void WriteEnsembleTest(bool readCalibration = false);
    void EvalEnsembleTest(bool writeCalibration = false);
    
    void QuickCalibration();
    void QuickSystematics();
    
    TH2F* Measure(Analysis* a);
};
