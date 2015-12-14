#ifndef TEMPLATE_DERIVATION_H
#define TEMPLATE_DERIVATION_H

#include "TTree.h"
#include "TChain.h"
#include "TH2F.h"

#include <string>

class RooWorkspace;
class RooDataSet;
class RooArgSet;
class RooRealVar;
class RooFormulaVar;
class RooFitResult;

class TemplateDerivation {
 public:
  TemplateDerivation();
  ~TemplateDerivation();

 private:
  TString selection_, samplePath_, fVar1_, fVar2_, fVar3_, fWeight_, fChannel_;
  std::string activeBranches_;

  const int maxPermutations_;
  const double maxMtop_;

  std::vector<double> jesValues_;
  std::vector<double> massValues_;

  enum channelID {
    kAllJets,
    kMuonJets,
    kElectronJets,
    kLeptonJets,
    kMaxChannels
  };
  int channelID_;
  RooWorkspace* workspace_;

  // do all the calculation
  void rooFitTopMass();

  // because RooFit only likes plain trees with standard data types (int, float,
  // double, ...)
  // the original tree has to be adapted for the new content
  TTree* modifiedTree(TChain* tree, int minComboType = -10,
                      int maxComboType = 10, bool isData = false);
  void UnknownChannelAbort();

  void addTemplateFunction(int templType, int comboType, int nComboTypes,
                           RooRealVar* mTop, RooRealVar* JES);
  RooFitResult* fitTemplate(int templType, int comboType, int nComboTypes);
  void plotResult(int templType, int comboType, int nComboTypes);
  void printResult(int templType, int comboType, int nComboTypes);
  void fillAlpha(std::vector<RooFormulaVar*>& alpha, int& h, RooArgSet argSet,
                 std::string add = "");

 public:
  std::vector<RooDataSet*> createDataSets(const RooArgSet& varSet);
  std::string constructFileName(double mass, double jes);
  static std::string constructTemplateName(double mass, double jes);
  unsigned int numVariables(TString startName) const;
  static TString constructName(const std::string& name, int i);
  static TString constructName(const std::string& name, int i, int j);
};

#endif /* TEMPLATE_DERIVATION_H */
