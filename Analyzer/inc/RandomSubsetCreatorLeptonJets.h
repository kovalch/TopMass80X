#ifndef RANDOMSUBSETCREATORLEPTONJETS_H_
#define RANDOMSUBSETCREATORLEPTONJETS_H_

#include "RandomSubsetCreator.h"

#include "TClonesArray.h"
#include "TFile.h"
#include "TH2F.h"
#include "TString.h"
#include "TTree.h"

class RandomSubsetCreatorLeptonJets : public RandomSubsetCreator {
public:
  RandomSubsetCreatorLeptonJets();
  virtual ~RandomSubsetCreatorLeptonJets();

private:
  TString selection_;
  TString samplePath_;
  TString fIdentifier_;
  TString fChannel_;
  double fLumi, fSig, fBDisc;
  std::string fWeight;

  TString fFile_;

  TFile* tmpFile_;

  int target, run, luminosityBlock, event, combi, nVertex, leptonId;
  double hadTopMass, hadWRawMass, leptonPt, leptonC, hitFitProb, deltaThetaHadWHadB, deltaThetaHadQHadQBar, PUWeight, PUWeightUp, PUWeightDown, muWeight, bWeight, bWeight_bTagSFUp, bWeight_bTagSFDown, bWeight_misTagSFUp, bWeight_misTagSFDown, MCWeight, mcWeight, hadQBCSV, hadQBarBCSV, hadBBCSV, lepBBCSV;
  double jetsPt[4];
  double pdfWeights[44];

  //TChain* fChain;

  TTree* fTree;
  TTree* fTreeTTmu;
  TTree* fTreeTTe;
  TTree* fTreeWmu;
  TTree* fTreeWe;
  TTree* fTreeSTmu;
  TTree* fTreeSTe;

  TTree* CreateRandomSubset();
  void DrawEvents(TTree* tempTree, double nEventsPE);
  TTree* PrepareTree(TString file);
};

#endif /* RANDOMSUBSETCREATORLEPTONJETS_H_ */
