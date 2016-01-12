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

  std::vector<double> jsfValues_;
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

  // because RooFit only likes plain trees with standard data types (int, float,
  // double, ...)
  // the original tree has to be adapted for the new content
  TTree* modifiedTree(TChain* tree, int minComboType = -10,
                      int maxComboType = 10, bool isData = false);
  void UnknownChannelAbort();
  void addTemplateFunction(const std::string& varType,
                           const std::string& comboType, RooRealVar* mTop,
                           RooRealVar* JES);
  RooFitResult* fitTemplate(const std::string& varType,
                            const std::string& comboType);
  void plotResult(const std::string& varType, const std::string& comboType);
  void printResult(const std::string& varType, const std::string& comboType);
  RooFormulaVar* createAlpha(const std::string& name, const RooArgSet& argSet,
                             const std::string& addTerm = "");
  std::vector<RooDataSet*> createDataSets(const RooArgSet* varSet);
  std::string constructFileName(double mass, double jes);
  unsigned int numVariables(const std::string& startName) const;
  void addTopCP();

 public:
  // do all the calculation
  void run();
  static std::string templateName(double mass, double jes);
  static std::string addInt(const std::string& name, int i);
  static std::string catName(const std::string& name, const std::string& cat1,
                          const std::string& cat2 = "",
                          const std::string& cat3 = "");
};

#endif /* TEMPLATE_DERIVATION_H */
