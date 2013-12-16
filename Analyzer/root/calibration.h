#ifndef ROOFITTOPMASS_C
#define ROOFITTOPMASS_C

#include "TTree.h"
#include "TChain.h"
#include "TH2F.h"

#include "RooFormulaVar.h"

#include <string>

class TopMassCalibration {
public:
  TopMassCalibration();
  ~TopMassCalibration();

private:
  TH2F *bTagEff_, *cTagEff_, *lTagEff_;

  TString selection_, samplePath_, fVar1_, fVar2_, fVar3_, fWeight_, fChannel_;
  std::string activeBranches_;

  int maxPermutations_;

  enum channelID {kAllJets, kMuonJets, kElectronJets, kLeptonJets, kMaxChannels};
  int channelID_;

  //int minComboType_ = 1;
  //int maxComboType_ = 1;

  //int templType_ = 0;

  bool doCalibration_;
  bool fitBackground_;
  //bool doMeasurement_;

  // do all the calculation
  void rooFitTopMass_();

  //// return the PU weights for the different samples
  //enum enumForPUWeights {kSummer11, kSummer11Plus05, kSummer11Minus05, kFall11, kFall11Plus05, kFall11Minus05, kFall10, kSummer12};
  //double calcPUWeight_(enum enumForPUWeights sample, short nPU);
  //
  //// calculate the probability of b-tagging one event with 2 b-tags
  //double eventBTagProbability_(std::vector<double> &oneMinusBEffies, std::vector<double> &oneMinusBMistags);
  //double calcBTagWeight_(int Njet, float * bTag, short * pdgId, TClonesArray * jets);

  // because RooFit only likes plain trees with standard data types (int, float, double, ...)
  // the original tree has to be adapted for the new content
  TTree* modifiedTree_(TChain *tree);
  TTree* modifiedTree_(TChain *tree, int comboType);
  TTree* modifiedTree_(TChain *tree, int minComboType, int maxComboType);
  TTree* modifiedTree_(TChain *tree, int minComboType, int maxComboType, bool isData);
  void UnknownChannelAbort();

  void fillAlpha(std::vector<RooFormulaVar*>& alpha, int& h, RooArgSet argSet);
};

#endif /* ROOFITTOPMASS_C */
