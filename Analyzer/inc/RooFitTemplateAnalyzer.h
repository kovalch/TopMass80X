#ifndef ROOFITTEMPLATEANALYZER_H
#define ROOFITTEMPLATEANALYZER_H

#include "MassAnalyzer.h"

#include "RooWorkspace.h"

#include "TString.h"

class RooFitTemplateAnalyzer : public MassAnalyzer { 
 public:
  RooFitTemplateAnalyzer(const std::string& identifier, TTree* tree);
  ~RooFitTemplateAnalyzer();
  void Analyze(const std::string& cuts, int i, int j);
    
 private:
  double QBTagProbability(double bDiscriminator);
  double mWnVertex();
  double mTnVertex();
  void Scan(const std::string& cuts, int i, int j, TString variables);
  
  static RooWorkspace* workspace;
};

#endif
