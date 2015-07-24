#include "FWCore/Utilities/interface/EDMException.h"
#include "TopAnalysis/TopUtils/plugins/myRecoMETProducer.h"

#include "DataFormats/METReco/interface/MET.h"
//#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

MyRecoMETProducer::MyRecoMETProducer(const edm::ParameterSet& cfg):
  mETSrc_( cfg.getParameter<edm::InputTag>("metSrc") )
{
  produces<std::vector<reco::MET> >();
}

MyRecoMETProducer::~MyRecoMETProducer()
{
}

void
MyRecoMETProducer::produce(edm::Event& evt, const edm::EventSetup& setup)
{
  edm::Handle<edm::View<reco::Candidate> > genMETs;
  evt.getByLabel(mETSrc_, genMETs);

  std::auto_ptr<std::vector<reco::MET> > vMETs(new std::vector<reco::MET>());

  for(unsigned int i = 0; i < genMETs->size(); ++i){
    vMETs->push_back(reco::MET(genMETs->at(i).p4(), reco::LeafCandidate::Point( 0, 0, 0 )));
  }
  // put our produced stuff in the event
  evt.put(vMETs);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE( MyRecoMETProducer );
