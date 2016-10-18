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
  virtual void endRun(const edm::Run&, const edm::EventSetup&);
  virtual void endJob();

  edm::Service<TreeRegistryService> trs;

  double mcWeight_;
  double topPtSF_[3];
  //edm::InputTag puSrc_;
  edm::EDGetTokenT<edm::View<PileupSummaryInfo> > puSrcToken_;
/*  edm::InputTag vertexSrc_;*/
  edm::EDGetTokenT<std::vector<reco::Vertex> > vertexSrc_;

  edm::InputTag puWeightSrc_;
  edm::InputTag puWeightUpSrc_;
  edm::InputTag puWeightDownSrc_;

  edm::InputTag jets_;
  edm::InputTag bWeightSrc_;
  edm::InputTag bWeightSrc_bTagSFUp_;
  edm::InputTag bWeightSrc_bTagSFDown_;
  edm::InputTag bWeightSrc_bTagCjetSFUp_;
  edm::InputTag bWeightSrc_bTagCjetSFDown_; 
  edm::InputTag bWeightSrc_misTagSFUp_;
  edm::InputTag bWeightSrc_misTagSFDown_;
  
  edm::InputTag bJESSrc_fNuUp_;
  edm::InputTag bJESSrc_fNuDown_;
  edm::InputTag bJESSrc_frag_;
  edm::InputTag bJESSrc_fragHard_;
  edm::InputTag bJESSrc_fragSoft_;

  edm::InputTag lepIDWeightSrc_; 
  edm::InputTag lepIDWeightSrcUp_;
  edm::InputTag lepIDWeightSrcDown_;
  
  edm::InputTag isoWeightSrc_;
  edm::InputTag isoWeightSrcUp_;
  edm::InputTag isoWeightSrcDown_;
  
  edm::InputTag triggerWeightSrc_;
  edm::InputTag triggerWeightSrcUp_;
  edm::InputTag triggerWeightSrcDown_;
  
  edm::InputTag muWeightSrc_;
  edm::InputTag elWeightSrc_;

  edm::InputTag genEventSrc_;
  edm::InputTag lheEventSrc_;
//   edm::InputTag ttEvent_;
  edm::InputTag ttInputTag;
  edm::EDGetTokenT<edm::View<TtEvent> > ttEvent_;
  bool savePDFWeights_;
  bool brCorrection_;
  bool showLHEweightTypes_;
  bool showLHEweightTypes2_;
  // THE WeightEvent to store the information
  WeightEvent* weight;
};

#endif
