#ifndef EventHypothesisAnalyzer_h
#define EventHypothesisAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "AnalysisDataFormats/TopObjects/interface/TtSemiLeptonicEvent.h"
#include "AnalysisDataFormats/TopObjects/interface/TtFullHadronicEvent.h"
#include "TopMass/TopEventTree/interface/TopEvent.h"


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
  
  // get vector of partons from the hypothesis
  std::vector<TLorentzVector> getPartons(const TtSemiLeptonicEvent *ttEvent, TtEvent::HypoClassKey hypoClassKey, unsigned int h);
  std::vector<TLorentzVector> getPartons(const TtFullHadronicEvent *ttEvent, TtEvent::HypoClassKey hypoClassKey, unsigned int h);
  
  // fills vector of genPartons, returns ID for the decay channel
  int fillGenPartons(const TtGenEvent *genEvent);

  edm::Service<TreeRegistryService> trs;

  edm::InputTag ttEvent_;
  edm::InputTag hypoClassKey_;
  
  //edm::InputTag leps_;
  //edm::InputTag mets_;
  //
  //bool data_;

  // max possible number of jets in events
  //const int kJetMAX_;

  // max possible number of permutations per event
  const unsigned int kMAXCombo_;

  // THE TopEvent to store the information
  TopEvent* top;
};

#endif
