#ifndef IDEOGRAMANALYZER_H
#define IDEOGRAMANALYZER_H

#include "MassAnalyzer.h"

class IdeogramAnalyzer : public MassAnalyzer { 
 public:
  IdeogramAnalyzer(const TString& identifier, TTree* tree) : MassAnalyzer(identifier, tree) {};
  void Analyze(const TString& cuts, int i, int j) { std::cout << "Function Analyze(const TString& cuts, int i, int j) not defined for Class IdeogramAnalyzer!" << std::endl; };
  void Analyze(const TString& cuts, int i, int j, po::variables_map vm);
    
 private:
  double QBTagProbability(double bDiscriminator);
  double mWnVertex();
  double mTnVertex();
  void Scan(const TString& cuts, int i, int j, double firstBinMass, double lastBinMass,
	    double resolMass, double firstBinJes, double lastBinJes, double resolJes, po::variables_map vm, bool fit2D = true);
  
};

#endif
