#ifndef IDEOGRAMANALYZER_H
#define IDEOGRAMANALYZER_H

#include "MassAnalyzer.h"

class IdeogramAnalyzer : public MassAnalyzer { 
 public:
  IdeogramAnalyzer(const std::string& identifier, TTree* tree) : MassAnalyzer(identifier, tree) {};
  void Analyze(const std::string& cuts, int i, int j);
    
 private:
  double QBTagProbability(double bDiscriminator);
  double mWnVertex();
  double mTnVertex();
  void Scan(const std::string& cuts, int i, int j, double firstBinMass, double lastBinMass,
	        double resolMass, double firstBinJes, double lastBinJes, double resolJes, bool fit2D = true);
  
};

#endif
