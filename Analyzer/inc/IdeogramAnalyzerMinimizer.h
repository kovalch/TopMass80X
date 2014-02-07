#ifndef IDEOGRAMANALYZERMINIMIZER_H
#define IDEOGRAMANALYZERMINIMIZER_H

#include "MassAnalyzer.h"

#include "DataSample.h"

class TF2;
class IdeogramCombLikelihood;

class IdeogramAnalyzerMinimizer : public MassAnalyzer {
 public:
  IdeogramAnalyzerMinimizer(const std::string& identifier, TTree* tree);
  ~IdeogramAnalyzerMinimizer();
  void Analyze(const std::string& cuts, int iBin, int jBin);
  void SetDataSample(const DataSample& sample) {sample_ = sample;}
    
 private:
  void Scan(const std::string& cuts, int iBin, int jBin);
  void NumericalMinimization();
  
  DataSample& sample_;
  //IdeogramCombLikelihood* fptr_;
  //TF2* combLikelihood_;
  int channelID_;
  double isFastSim_;
  double shapeSystematic_;
  double permutationFractionSystematic_;
  //std::string topBranchName_;
  
  std::vector<std::vector<IdeogramCombLikelihood*>> eventFunctions_;

};

#endif
