#ifndef TOPMASS_H
#define TOPMASS_H

#include <vector>
#include <boost/program_options.hpp>

#include "TString.h"

#include "Analysis.h"

namespace po = boost::program_options;

/*
struct massPoint {
  double genMass;
  double genJES;
  double genLumi;
  TString identifier;
  TString fileName;
  Analysis* analysis;
  
  TH1D* hMass;
  TH1D* hMassError;
  TH1D* hMassPull;
  TH2F* h2Mass;
  TH2F* h2MassError;
  TH3F* h3Mass;
  TH3F* h3MassError;
  TH3F* h3MassPull;
  TH3F* h3JES;
  TH3F* h3JESError;
  TH3F* h3JESPull;
};
*/

class TopMass {
  private:
    TString fMethod_;
    int fBins_;
    double fLumi_;
    
    std::vector< std::vector<Analysis*> > calibrationAnalyses;
    
    std::vector<Analysis*> analyses;
    std::vector<Analysis*>::iterator iAnalysis;
    
    Analysis* aSim;
    
    double fCalibFitParameter[6][6][2];
    double fCalibFitParError[6][6][2];
    
    //void LoadXML();
    bool fexists(const char *filename);
  
  public:
    TopMass(po::variables_map vm);
    
    void WriteEnsembleTest(po::variables_map vm);
    
    void QuickCalibration(po::variables_map vm);
    void QuickSystematics(po::variables_map vm);

    //void AvoidCompilerWarning(){}
};

#endif
