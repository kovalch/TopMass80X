#ifndef IDEOGRAMANALYZERMINIMIZER_H
#define IDEOGRAMANALYZERMINIMIZER_H

#include "MassAnalyzer.h"

#include "DataSample.h"

class TF2;
class IdeogramCombLikelihood;
namespace ROOT{
  namespace Math{
  class Minimizer;
  }
}

class IdeogramAnalyzerMinimizer : public MassAnalyzer {
 public:
  IdeogramAnalyzerMinimizer(const std::string& identifier, TTree* tree);
  ~IdeogramAnalyzerMinimizer();
  void Analyze(const std::string& cuts, int iBin, int jBin);
  void SetDataSample(const DataSample& sample) {sample_ = sample;}
    
 private:
  void Scan(const std::string& cuts, int iBin, int jBin);
  void NumericalMinimization();
  enum allowedVariables{kMass, kJES, kFSig, kFCP};
  void IterateVariableCombinations(ROOT::Math::Minimizer* min, std::vector<allowedVariables> toFit, unsigned int start = 0);
  void PlotResult(ROOT::Math::Minimizer* min, allowedVariables x = kMass, allowedVariables y = kJES);

  void CleanUp();
  
  DataSample& sample_;
  int channelID_;
  double isFastSim_;
  double shapeSystematic_;
  double permutationFractionSystematic_;
  //std::string topBranchName_;
  
  std::vector<std::vector<IdeogramCombLikelihood*>> eventFunctions_;

};

#endif
