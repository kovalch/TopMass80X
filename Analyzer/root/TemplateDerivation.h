#ifndef TEMPLATE_DERIVATION_H
#define TEMPLATE_DERIVATION_H

#include "TTree.h"
#include "TChain.h"
#include "TH2F.h"

#include "RooFormulaVar.h"

#include <string>

class TemplateDerivation {
 public:
  TemplateDerivation();
  ~TemplateDerivation();

 private:
  TString selection_, samplePath_, fVar1_, fVar2_, fVar3_, fWeight_, fChannel_;
  std::string activeBranches_;

  const int maxPermutations_;
  const double maxMtop_;

  enum channelID {
    kAllJets,
    kMuonJets,
    kElectronJets,
    kLeptonJets,
    kMaxChannels
  };
  int channelID_;

  bool doCalibration_;
  bool fitBackground_;

  // do all the calculation
  void rooFitTopMass_();

  // because RooFit only likes plain trees with standard data types (int, float,
  // double, ...)
  // the original tree has to be adapted for the new content
  TTree* modifiedTree_(TChain* tree, int comboType);
  TTree* modifiedTree_(TChain* tree, int minComboType = -10,
                       int maxComboType = 10, bool isData = false);
  void UnknownChannelAbort();

  void fillAlpha(std::vector<RooFormulaVar*>& alpha, int& h, RooArgSet argSet,
                 std::string add = "");
  std::vector<RooDataSet*> createDataSets(std::vector<double> masses,
                                          std::vector<double> JESes,
                                          const RooArgSet& varSet);
  std::string constructFileName(double mass, double jes);
};

#endif /* TEMPLATE_DERIVATION_H */
