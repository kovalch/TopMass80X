#ifndef IDEOGRAMANALYZERNEWINTERFACE_H
#define IDEOGRAMANALYZERNEWINTERFACE_H

#include "MassAnalyzer.h"

#include "IdeogramCombLikelihood.h"

#include "TF2.h"

class IdeogramAnalyzerNewInterface : public MassAnalyzer {
 public:
  IdeogramAnalyzerNewInterface(const TString& identifier, TTree* tree);
  ~IdeogramAnalyzerNewInterface();
  void Analyze(const TString& cuts, int i, int j);
    
 private:
  double QBTagProbability(double bDiscriminator);
  double mWnVertex();
  double mTnVertex();
  void Scan(const TString& cuts, int i, int j, double firstBinMass, double lastBinMass,
	        double resolMass, double firstBinJes, double lastBinJes, double resolJes, bool fit2D = true);
  
  IdeogramCombLikelihood* fptr_;
  TF2* combLikelihood_;
  int channelID_;
};

#endif
