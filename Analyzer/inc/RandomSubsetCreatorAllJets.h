/*
 * RandomSubsetCreatorAllJets.h
 *
 *  Created on: Oct 23, 2012
 *      Author: eschliec
 */

#ifndef RANDOMSUBSETCREATORALLJETS_H_
#define RANDOMSUBSETCREATORALLJETS_H_

#include "RandomSubsetCreator.h"

//#include "boost/variant.hpp"
#include "boost/program_options.hpp"

#include "TH2F.h"
#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TClonesArray.h"

namespace po = boost::program_options;

class RandomSubsetCreatorAllJets : public RandomSubsetCreator {
public:
  RandomSubsetCreatorAllJets(po::variables_map vm);
  virtual ~RandomSubsetCreatorAllJets();

private:
  TString selection_;
  TString samplePath_;

  double fLumi_;

  TString fIdentifier_;
  TString fFile_;

  TFile* tmpFile_;

  TH2F* bTagEff_;
  TH2F* cTagEff_;
  TH2F* lTagEff_;

  //typedef boost::variant< short, int, float, double, unsigned int, unsigned short*, short*, double*, TClonesArray*> variableTypes;

  //std::map<TString, std::map<TString, variableTypes> > _variables;

  //static xml::XMLDocument* _config;

  void FetchBTagEfficiencyHistograms();
  //void ReadConfigFromXMLFile();

  TTree* CreateRandomSubset();

  //variableTypes& GetVariable(TString variableName);

  //void SetBranchStatuses(TTree* tree);
  //void SetBranchAddresses(TTree* tree);

  // return the PU weights for the different samples
  enum enumForPUWeights {kSummer11, kSummer11Plus05, kSummer11Minus05, kFall11, kFall11Plus05, kFall11Minus05, kFall10};
  double calcPUWeight_(enum enumForPUWeights sample, short nPU);
  double calcPDFWeight_(int whichPDFUncertainty, bool upVariation, double x1, int id1, float Q, double x2, int id2);
  double eventBTagProbability_(std::vector<double> &oneMinusBEffies, std::vector<double> &oneMinusBMistags, bool verbose = false);
  double calcBTagWeight_(int Njet=0, short* pdgId=0, TClonesArray* jets=0);

  void AddWeights(TTree* tempTree, bool isData=false);
};

#endif /* RANDOMSUBSETCREATORALLJETS_H_ */
