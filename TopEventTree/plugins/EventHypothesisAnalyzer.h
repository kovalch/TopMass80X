#ifndef EventHypothesisAnalyzer_h
#define EventHypothesisAnalyzer_h

#include "TTree.h"

#include <vector>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "AnalysisDataFormats/TopObjects/interface/TtSemiLeptonicEvent.h"
#include "AnalysisDataFormats/TopObjects/interface/TtFullHadronicEvent.h"
#include "TopMass/TopEventTree/interface/TopEvent.h"

//typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >   PxPyPzEVector;
//typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double> > PtEtaPhiEVector;

class EventHypothesisAnalyzer : public edm::EDAnalyzer {

 public:

  explicit EventHypothesisAnalyzer(const edm::ParameterSet&);
  ~EventHypothesisAnalyzer();
  
  enum semiLepParticles{ TTBar, HadTop, LepTop, HadW, LepW, HadB, LightQ, LightQBar, LepB, Lepton, Neutrino };
  enum fullHadParticles{ TTBar0, Top1, Top2, W1, W2, B1, LightQ1, LightQBar1, B2, LightQ2, LightQBar2 };
  
 private:

  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  
  std::vector<TLorentzVector> getPartons(const TtSemiLeptonicEvent *ttEvent, TtEvent::HypoClassKey hypoClassKey, unsigned int h);
  std::vector<TLorentzVector> getPartons(const TtFullHadronicEvent *ttEvent, TtEvent::HypoClassKey hypoClassKey, unsigned int h);
  
  edm::Service<TreeRegistryService> trs;

  TTree* eventTree;
  TBranch* hypoBranch;
  
  edm::InputTag ttEvent_;
  edm::InputTag hypoClassKey_;
  
  edm::InputTag jets_;
  edm::InputTag allJets_;
  edm::InputTag noPtEtaJets_;
  edm::InputTag leps_;
  edm::InputTag mets_;

  edm::InputTag PUSrc_;
  edm::InputTag VertexSrc_;

  edm::InputTag PUWeightSrc_;
  edm::InputTag PUWeightUpSrc_;
  edm::InputTag PUWeightDownSrc_;

  edm::InputTag PUAWeightSrc_;
  edm::InputTag PUAWeightUpSrc_;
  edm::InputTag PUAWeightDownSrc_;

  edm::InputTag PUBWeightSrc_;
  edm::InputTag PUBWeightUpSrc_;
  edm::InputTag PUBWeightDownSrc_;

  edm::InputTag PUABWeightSrc_;
  edm::InputTag PUABWeightUpSrc_;
  edm::InputTag PUABWeightDownSrc_;

  edm::InputTag bWeightSrc_;
  edm::InputTag bWeightSrc_bTagSFUp_;
  edm::InputTag bWeightSrc_bTagSFDown_;
  edm::InputTag bWeightSrc_misTagSFUp_;
  edm::InputTag bWeightSrc_misTagSFDown_;
  
  edm::InputTag muWeightSrc_;
  edm::InputTag mcWeightSrc_;
  
  bool savePDFWeights_;
  bool data_;

  // max possible number of jets in events
  const int kJetMAX_;

  // max possible number of permutations per event
  const unsigned int kMAXCombo_;

  int nJet;
  int permutation[4];
  int noPtEtaJetMultiplicity;
  int bottomSSVJetMultiplicity;
  int bottomCSVJetMultiplicity;
  
  
  TopEvent* top;
};

#endif
