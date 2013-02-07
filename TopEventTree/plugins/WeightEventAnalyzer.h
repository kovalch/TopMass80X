#ifndef WeightEventAnalyzer_h
#define WeightEventAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "TopMass/TopEventTree/interface/WeightEvent.h"


class WeightEventAnalyzer : public edm::EDAnalyzer {

 public:

  explicit WeightEventAnalyzer(const edm::ParameterSet&);
  ~WeightEventAnalyzer();

 private:

  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  edm::Service<TreeRegistryService> trs;

  edm::InputTag mcWeightSrc_;

  edm::InputTag puSrc_;
  edm::InputTag vertexSrc_;

  edm::InputTag puWeightSrc_;
  edm::InputTag puWeightUpSrc_;
  edm::InputTag puWeightDownSrc_;

  edm::InputTag jets_;
  edm::InputTag bWeightSrc_;
  edm::InputTag bWeightSrc_bTagSFUp_;
  edm::InputTag bWeightSrc_bTagSFDown_;
  edm::InputTag bWeightSrc_misTagSFUp_;
  edm::InputTag bWeightSrc_misTagSFDown_;

  edm::InputTag muWeightSrc_;
  edm::InputTag elWeightSrc_;

  edm::InputTag genEventSrc_;
  bool savePDFWeights_;

  // THE WeightEvent to store the information
  WeightEvent* weight;
};

#endif
