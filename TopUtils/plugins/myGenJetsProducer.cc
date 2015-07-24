#include "FWCore/Utilities/interface/EDMException.h"
#include "TopAnalysis/TopUtils/plugins/myGenJetsProducer.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/JetReco/interface/GenJet.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

MyGenJetsProducer::MyGenJetsProducer(const edm::ParameterSet& cfg):
  jetsSrc_( cfg.getParameter<edm::InputTag>("genJets") ),
  leptonSrc_( cfg.getParameter<edm::InputTag>("genParticles") )
{
  produces<std::vector<reco::GenJet> >();
}

MyGenJetsProducer::~MyGenJetsProducer()
{
}

void
MyGenJetsProducer::produce(edm::Event& evt, const edm::EventSetup& setup)
{
  edm::Handle<std::vector<reco::GenJet> > genJets;
  evt.getByLabel(jetsSrc_, genJets);

  edm::Handle<edm::View<reco::Candidate> > genParticles;
  evt.getByLabel(leptonSrc_, genParticles);

  std::auto_ptr<std::vector<reco::GenJet> > vGenJets(new std::vector<reco::GenJet>());

  for(unsigned int i = 0; i < genJets->size(); ++i){
    bool matched=false;
    for(unsigned int j = 0; j < genParticles->size(); ++j){
        //std::cout << genParticles->at(j).pdgId() << " " << genParticles->at(j).pt() << " " << genParticles->at(j).eta() << " " << genParticles->at(j).phi() << std::endl;
        if(deltaR2(genParticles->at(j),genJets->at(i)) < 0.5*0.5){
            matched=true;
            break;
        }
    }
    if (!matched){
        vGenJets->push_back(reco::GenJet(genJets->at(i).p4(), genJets->at(i).getSpecific(), genJets->at(i).getJetConstituents()));
    }
  }
  // put our produced stuff in the event
  evt.put(vGenJets);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE( MyGenJetsProducer );
