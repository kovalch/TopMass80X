#ifndef ROOFITTEMPLATEANALYZER_H
#define ROOFITTEMPLATEANALYZER_H

#include "MassAnalyzer.h"

#include "RooWorkspace.h"

class RooFitTemplateAnalyzer : public MassAnalyzer { 
 public:
  RooFitTemplateAnalyzer(const TString& identifier, TTree* tree);
  ~RooFitTemplateAnalyzer();
  void Analyze(const TString& cuts, int i, int j) { std::cout << "Function Analyze(const TString& cuts, int i, int j) not defined for Class RooFitTemplateAnalyzer!" << std::endl; };
  void Analyze(const TString& cuts, int i, int j, po::variables_map vm);
    
 private:
  double QBTagProbability(double bDiscriminator);
  double mWnVertex();
  double mTnVertex();
  void Scan(const TString& cuts, int i, int j, po::variables_map vm, TString variables);
  
  static RooWorkspace* workspace;
};

#endif
