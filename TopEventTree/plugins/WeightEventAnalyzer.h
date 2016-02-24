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

  double mcWeight_;

  edm::InputTag puSrc_;
/*  edm::InputTag vertexSrc_;*/
  edm::EDGetTokenT<std::vector<reco::Vertex> > vertexSrc_;

  edm::InputTag puWeightSrc_;
  edm::InputTag puWeightUpSrc_;
  edm::InputTag puWeightDownSrc_;

  edm::InputTag jets_;
  edm::InputTag bWeightSrc_;
  edm::InputTag bWeightSrc_bTagSFUp_;
  edm::InputTag bWeightSrc_bTagSFDown_;
  edm::InputTag bWeightSrc_misTagSFUp_;
  edm::InputTag bWeightSrc_misTagSFDown_;
  
  edm::InputTag bJESSrc_fNuUp_;
  edm::InputTag bJESSrc_fNuDown_;
  edm::InputTag bJESSrc_frag_;
  edm::InputTag bJESSrc_fragHard_;
  edm::InputTag bJESSrc_fragSoft_;

  edm::InputTag triggerWeightSrc_;

  edm::InputTag muWeightSrc_;
  edm::InputTag elWeightSrc_;

  edm::InputTag genEventSrc_;
  edm::InputTag lheEventSrc_;
//   edm::InputTag ttEvent_;
  edm::InputTag ttInputTag;
  edm::EDGetTokenT<edm::View<TtEvent> > ttEvent_;
  bool savePDFWeights_;
  bool brCorrection_;

  // THE WeightEvent to store the information
  WeightEvent* weight;
};

#endif
