#include "FWCore/Utilities/interface/EDMException.h"
#include "TopAnalysis/TopUtils/plugins/myGenElectronToMuonConverter.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

MyGenElectronToMuonConverter::MyGenElectronToMuonConverter(const edm::ParameterSet& cfg):
  electronSrc_( cfg.getParameter<edm::InputTag>("electronSrc") )
{
  produces<std::vector<reco::GenParticle> >();
}

MyGenElectronToMuonConverter::~MyGenElectronToMuonConverter()
{
}

void
MyGenElectronToMuonConverter::produce(edm::Event& evt, const edm::EventSetup& setup)
{
  edm::Handle<edm::View<reco::Candidate> > genElectrons;
  evt.getByLabel(electronSrc_, genElectrons);

  std::auto_ptr<std::vector<reco::GenParticle> > vElectrons(new std::vector<reco::GenParticle>());

  for(unsigned int i = 0; i < genElectrons->size(); ++i){
    vElectrons->push_back(reco::GenParticle(*dynamic_cast<const reco::LeafCandidate*>(&genElectrons->at(i))));
    if     (std::abs(genElectrons->at(i).pdgId()) == 11) vElectrons->back().setPdgId(-13*genElectrons->at(i).charge());
    else if(std::abs(genElectrons->at(i).pdgId()) == 13) vElectrons->back().setPdgId(-11*genElectrons->at(i).charge());
  }
  // put our produced stuff in the event
  evt.put(vElectrons);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE( MyGenElectronToMuonConverter );
