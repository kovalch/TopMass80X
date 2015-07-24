#include "FWCore/Utilities/interface/EDMException.h"
#include "TopAnalysis/TopUtils/plugins/myRecoMuonProducer.h"

#include "DataFormats/MuonReco/interface/Muon.h"
//#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

MyRecoMuonProducer::MyRecoMuonProducer(const edm::ParameterSet& cfg):
  muonSrc_( cfg.getParameter<edm::InputTag>("muonSrc") )
{
  produces<std::vector<reco::Muon> >();
}

MyRecoMuonProducer::~MyRecoMuonProducer()
{
}

void
MyRecoMuonProducer::produce(edm::Event& evt, const edm::EventSetup& setup)
{
  edm::Handle<edm::View<reco::Candidate> > genMuons;
  evt.getByLabel(muonSrc_, genMuons);

  std::auto_ptr<std::vector<reco::Muon> > vMuons(new std::vector<reco::Muon>());

  for(unsigned int i = 0; i < genMuons->size(); ++i){
    vMuons->push_back(reco::Muon(genMuons->at(i).charge(), genMuons->at(i).p4()));
    vMuons->back().setPdgId(genMuons->at(i).pdgId());
  }
  // put our produced stuff in the event
  evt.put(vMuons);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE( MyRecoMuonProducer );
