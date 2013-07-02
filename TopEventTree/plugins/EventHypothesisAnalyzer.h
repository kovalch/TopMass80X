#ifndef EventHypothesisAnalyzer_h
#define EventHypothesisAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "AnalysisDataFormats/TopObjects/interface/TtSemiLeptonicEvent.h"
#include "AnalysisDataFormats/TopObjects/interface/TtFullHadronicEvent.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "TopMass/TopEventTree/interface/TopEvent.h"


class EventHypothesisAnalyzer : public edm::EDAnalyzer {

 public:

  explicit EventHypothesisAnalyzer(const edm::ParameterSet&);
  ~EventHypothesisAnalyzer();
  
 private:

  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  
  // get vector of partons from the hypothesis
  std::vector<TLorentzVector> getPartons(const TtSemiLeptonicEvent *ttEvent, TtEvent::HypoClassKey hypoClassKey, unsigned int h);
  std::vector<TLorentzVector> getPartons(const TtFullHadronicEvent *ttEvent, TtEvent::HypoClassKey hypoClassKey, unsigned int h);
  // needed in fullHad channel for reco masses
  std::vector<TLorentzVector> getPartons(const std::vector<pat::Jet> *jets);

  // fills vector of genPartons, returns ID for the decay channel
  int fillGenPartons(const TtGenEvent *genEvent);

  short comboTypeFullHad();
  short comboTypeSemiLep();
  short comboTypeIDCalculatorFullHad();
  short comboTypeIDCalculatorSemiLep();
  short comboTypeAlgo(std::vector<int> jetIndexFit, std::vector<int> jetIndexGen);
  std::vector<short> comboTypeAlgoInverted(std::vector<int> jetIndexGen, short comboType);

  edm::Service<TreeRegistryService> trs;

  edm::InputTag ttEvent_;
  edm::InputTag hypoClassKey_;
  edm::InputTag ttEventGen2_;
  edm::InputTag jets_;
  
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
