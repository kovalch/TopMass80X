#ifndef RANDOMSUBSETCREATORNEWINTERFACE_H_
#define RANDOMSUBSETCREATORNEWINTERFACE_H_

#include "Helper.h"
#include "RandomSubsetCreator.h"
#include "DataSample.h"

#include <string>

class TRandom3;

class RandomSubsetCreatorNewInterface : public RandomSubsetCreator {
public:
  RandomSubsetCreatorNewInterface(const std::vector<float>& v);
  virtual ~RandomSubsetCreatorNewInterface();
  const DataSample& GetDataSample() const { return subset_; }

private:
  const std::string selection_;
  const std::string selectionLept_;
  const std::string selectionJets_;
  const std::string samplePath_;
  const std::string samplePathLept_;
  const std::string samplePathJets_;
  const std::string fIdentifier_;
  const std::string fVar1_;
  const std::string fVar2_;
  const std::string fVar3_;
  const std::string fVar4_;
  const std::string fWeight_;
  const std::string activeBranches_;
  const std::string fBinning_;
  const std::vector<float> vBinning_;
  const double fLumi_, fSig_; //, fBDisc_;
  double fLumiLept_, fLumiJets_;
  double fSigLept_, fSigJets_;
  const int maxPermutations_;
  int channelID_;

  TRandom3* random_;

  std::vector<DataSample> events_;
  DataSample mergedsample_;
  DataSample subset_;

//always return 0! get the output as DataSample over GetDataSample()
  TTree* CreateRandomSubset();
  void DrawEvents(const DataSample& sample, double nEventsPE);
  void PrepareEvents(const std::string& file, const Helper::ChannelID currentID = Helper::kMaxChannels, double sampleFactor = 1.);
};

#endif /* RANDOMSUBSETCREATORNEWINTERFACE_H_ */
