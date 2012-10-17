#include "MassAnalyzer.h"
#include "TF1.h"
#include "TCanvas.h"

#include <iostream>

class MVAAnalyzer : public MassAnalyzer {
  private:
    TF1* gaus;
  
  public:
    MVAAnalyzer(const TString& identifier, TTree* tree) : MassAnalyzer(identifier, tree), gaus(0) {};
    void Analyze(const TString& cuts, int i, int j);
    void Analyze(const TString& cuts, int i, int j, po::variables_map vm) { std::cout << "Function Analyze(const TString& cuts, int i, int j, po::variables_map vm) not defined for Class GenMatchAnalyzer!" << std::endl; };
    
    double GetMass();
};
