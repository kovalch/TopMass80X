#ifndef IDEOGRAMANALYZERNEWINTERFACE_H
#define IDEOGRAMANALYZERNEWINTERFACE_H

#include "MassAnalyzer.h"

#include "DataSample.h"

class TF2;
class IdeogramCombLikelihood;

class IdeogramAnalyzerNewInterface : public MassAnalyzer {
 public:
  IdeogramAnalyzerNewInterface(const std::string& identifier, TTree* tree);
  ~IdeogramAnalyzerNewInterface();
  void Analyze(const std::string& cuts, int i, int j);
  void SetDataSample(const DataSample& sample) {sample_ = sample;}
    
 private:
  double QBTagProbability(double bDiscriminator);
  double mWnVertex();
  double mTnVertex();
  void Scan(const std::string& cuts, int i, int j, double firstBinMass, double lastBinMass,
	        double resolMass, double firstBinJes, double lastBinJes, double resolJes, bool fit2D = true);
  
  DataSample& sample_;
  IdeogramCombLikelihood* fptr_;
  TF2* combLikelihood_;
  int channelID_;
  double pullWidth_;
  double isFastSim_;
  double shapeSystematic_;
  double permutationFractionSystematic_;
  std::string topBranchName_;

};

#endif
