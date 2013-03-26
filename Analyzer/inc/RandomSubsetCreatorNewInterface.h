#ifndef RANDOMSUBSETCREATORNEWINTERFACE_H_
#define RANDOMSUBSETCREATORNEWINTERFACE_H_

#include "RandomSubsetCreator.h"

#include "TClonesArray.h"
#include "TFile.h"
#include "TH2F.h"
#include "TRandom3.h"
#include "TString.h"
#include "TTree.h"

#include "TopMass/TopEventTree/interface/TopEvent.h"
#include "TopMass/TopEventTree/interface/WeightEvent.h"

class RandomSubsetCreatorNewInterface : public RandomSubsetCreator {
public:
  RandomSubsetCreatorNewInterface();
  virtual ~RandomSubsetCreatorNewInterface();

private:
  const TString selection_;
  const TString samplePath_;
  const TString fIdentifier_;
  //const TString fChannel_;
  const std::string fWeight_;
  const double fLumi_, fSig_, fBDisc_;
  int channelID_;

  TString fFile_;

  TopEvent* topEvent_;
  WeightEvent* weightEvent_;

  std::vector<TTree*> fTreesSig_;
  std::vector<TTree*> fTreesBkg_;
  TTree* fTree_;
  TRandom3* random_;

  TTree* CreateRandomSubset();
  void DrawEvents(TTree* tempTree, double nEventsPE);
  TTree* PrepareTree(TString file);
};

#endif /* RANDOMSUBSETCREATORNEWINTERFACE_H_ */
