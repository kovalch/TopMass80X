#ifndef JetEventMixer_h
#define JetEventMixer_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Sources/interface/VectorInputSource.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

class JetEventMixer : public edm::EDProducer {

public:
  explicit JetEventMixer(const edm::ParameterSet&);
  ~JetEventMixer();

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  void getEvents(edm::Event& evt);
  void processOneEvent(edm::EventPrincipal const& eventPrincipal, edm::Event& evt);
  void putOneEvent(edm::Event& evt);
  void cleanUp();

  edm::VectorInputSource *eventSrc_;

  const unsigned int nMix_;
  std::vector<unsigned int> combo_;

  std::vector<std::vector<pat::Jet         > > oriPatJets_;
  //std::vector<std::vector<reco::GenJet     > > oriGenJets_;
  std::vector<std::vector<reco::PFCandidate> > oriPFCandidates_;
  std::vector<std::vector<PileupSummaryInfo> > oriPUInfos_;
};

#endif
