#ifndef IDEOGRAMANALYZERNEWINTERFACE_H
#define IDEOGRAMANALYZERNEWINTERFACE_H

#include "MassAnalyzer.h"

class IdeogramAnalyzerNewInterface : public MassAnalyzer {
 public:
  IdeogramAnalyzerNewInterface(const TString& identifier, TTree* tree) : MassAnalyzer(identifier, tree) {};
  void Analyze(const TString& cuts, int i, int j);
    
 private:
  double QBTagProbability(double bDiscriminator);
  double mWnVertex();
  double mTnVertex();
  void Scan(const TString& cuts, int i, int j, double firstBinMass, double lastBinMass,
	        double resolMass, double firstBinJes, double lastBinJes, double resolJes, bool fit2D = true);
  
};

#endif
