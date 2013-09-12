#ifndef RANDOMSUBSETCREATORNEWINTERFACE_H_
#define RANDOMSUBSETCREATORNEWINTERFACE_H_

#include "RandomSubsetCreator.h"
#include "DataSample.h"

class TRandom3;

class RandomSubsetCreatorNewInterface : public RandomSubsetCreator {
public:
  RandomSubsetCreatorNewInterface();
  virtual ~RandomSubsetCreatorNewInterface();
  const DataSample& GetDataSample() const { return subset_; }

private:
  const std::string selection_;
  const std::string samplePath_;
  const std::string fIdentifier_;
  const std::string fVar1_;
  const std::string fVar2_;
  const std::string fVar3_;
  const std::string fWeight_;
  const std::string activeBranches_;
  const double fLumi_, fSig_, fBDisc_;
  int channelID_;

  TRandom3* random_;

  std::vector<DataSample> events_;
  DataSample subset_;

  TTree* CreateRandomSubset();
  void DrawEvents(const DataSample& sample, double nEventsPE);
  void PrepareEvents(const TString& file);
};

#endif /* RANDOMSUBSETCREATORNEWINTERFACE_H_ */
