#ifndef RANDOMSUBSETCREATORNEWINTERFACE_H_
#define RANDOMSUBSETCREATORNEWINTERFACE_H_

#include "RandomSubsetCreator.h"
#include "DataSample.h"

//#include "TChain.h"
//#include "TClonesArray.h"
//#include "TFile.h"
//#include "TH2F.h"
#include "TRandom3.h"
#include "TString.h"

//#include "TopMass/TopEventTree/interface/TopEvent.h"
//#include "TopMass/TopEventTree/interface/WeightEvent.h"

class RandomSubsetCreatorNewInterface : public RandomSubsetCreator {
public:
  RandomSubsetCreatorNewInterface();
  virtual ~RandomSubsetCreatorNewInterface();

private:
  const TString selection_;
  const TString samplePath_;
  const TString fIdentifier_;
  //const TString fChannel_;
  const std::string fVar1_;
  const std::string fVar2_;
  const std::string fVar3_;
  const std::string fWeight_;
  const double fLumi_, fSig_, fBDisc_;
  int channelID_;

  //TString fFileName_;

  //TopEvent* topEvent_;
  //WeightEvent* weightEvent_;

  //std::vector<TChain*> fChainsSig_;
  //std::vector<TChain*> fChainsBkg_;
  //TChain* fChain_;
  TRandom3* random_;

  std::vector<DataSample> events_;
  DataSample subset_;

  TTree* CreateRandomSubset();
  DataSample GetDataSample() { return subset_; }
  void DrawEvents(const DataSample& sample, double nEventsPE);
  void PrepareEvents(TString file);
};

#endif /* RANDOMSUBSETCREATORNEWINTERFACE_H_ */
