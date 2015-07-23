#include "FWCore/Utilities/interface/EDMException.h"
#include "TopAnalysis/TopUtils/plugins/MyBTagProducer.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

MyBTagProducer::MyBTagProducer(const edm::ParameterSet& cfg):
  jetSrc_    ( cfg.getParameter<edm::InputTag>("jetSrc") ),
  bDiscValue_( cfg.getParameter<double>("bDiscValue") ),
  jetsToTag_ ( cfg.getParameter<std::vector<int> >("jetsToTag") )
{
  produces<reco::JetTagCollection>();
}

MyBTagProducer::~MyBTagProducer()
{
}

void
MyBTagProducer::produce(edm::Event& evt, const edm::EventSetup& setup)
{
  edm::Handle<edm::View<reco::Jet> > jetSrc;
  evt.getByLabel(jetSrc_, jetSrc);

  if(jetSrc.operator->()->size() == 0) return;
  edm::RefToBase<reco::Jet> jetRefToBase(jetSrc.operator->()->refAt(0));
  edm::RefToBaseProd<reco::Jet> jetRefToBaseProd(jetRefToBase);
  std::auto_ptr<reco::JetTagCollection> jetTagCollection(new reco::JetTagCollection(jetRefToBaseProd));

  for(int iJet = 0, nJets = (int)jetSrc->size(); iJet < nJets; ++iJet){
    const reco::Jet* jet = &jetSrc->at(iJet);
    bool isBJet = false;
    if(jetsToTag_.size() == 0){ 
      for(int i = 0, l = jet->numberOfDaughters(); i < l; ++i){
	const reco::Candidate* dau = containsBHadron(jet->daughter(i));
	if(dau){
	  if(deltaR(dau->eta(), dau->phi(), jet->eta(), jet->phi()) < 0.3){
	    //std::cout << "!!! Tagged !!! JET: "  << jet->pt() << " " << jet->eta() << " " << jet->phi() << std::endl;
	    isBJet = true;
	    break;
	  }
	}
      }
    }
    else{
      for(std::vector<int>::const_iterator it = jetsToTag_.begin(); it != jetsToTag_.end(); ++it){
	if(iJet == *it) {
	  isBJet = true;
	  break;
	}
      }
    }
    if(isBJet) (*jetTagCollection)[jetSrc.operator->()->refAt(iJet)] = bDiscValue_;
    else       (*jetTagCollection)[jetSrc.operator->()->refAt(iJet)] = 0.0;
  }

  evt.put(jetTagCollection);
}

const reco::Candidate* //bool
MyBTagProducer::containsBHadron(const reco::Candidate* consti)
{
  //std::cout << consti->pdgId() << " ";
  if(std::abs(consti->pdgId()        / 1000) == 5 ||
     std::abs(consti->pdgId() % 1000 /  100) == 5 ){
    //std::cout << " found b-hadron! " << consti->pt() << " " << consti->eta() << " " << consti->phi() << std::endl;
    return consti; //true;
  }
  else{
    for(int i = 0, l = consti->numberOfMothers(); i < l; ++i){
      const reco::Candidate* dau = consti->mother(i);
      return containsBHadron(dau);
      //if(containsBHadron(dau))
      //	return true;
      //else
      //	return false;
    }
  }
  return 0; //false;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE( MyBTagProducer );
