// -*- C++ -*-
//
// Package:    HiggsGenEventReco
// Class:      HiggsGenEventReco
// 
/**\class HiggsGenEventReco HiggsGenEventReco.cc TopAnalysis/HiggsUtils/plugins/HiggsGenEventReco.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Johannes Hauk,,,DESY
//         Created:  Mon Jan 14 18:52:07 CET 2013
// $Id: HiggsGenEventReco.cc,v 1.2 2013/01/28 16:39:33 hauk Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "TopAnalysis/HiggsUtils/interface/HiggsGenEvent.h"





//
// class declaration
//

class HiggsGenEventReco : public edm::EDProducer {
   public:
      explicit HiggsGenEventReco(const edm::ParameterSet&);
      ~HiggsGenEventReco();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
      
      edm::InputTag decayChainSource_;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
HiggsGenEventReco::HiggsGenEventReco(const edm::ParameterSet& iConfig):
decayChainSource_(iConfig.getParameter<edm::InputTag>("src"))
{
  produces<HiggsGenEvent>();
}


HiggsGenEventReco::~HiggsGenEventReco()
{
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
HiggsGenEventReco::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::GenParticleCollection> parts;
  iEvent.getByLabel(decayChainSource_, parts);
  
  // add Higgs decay tree
  reco::GenParticleRefProd cands( parts );
  
  //add genEvtHiggs to the output stream
  HiggsGenEvent* genEvt = new HiggsGenEvent( cands );
  std::auto_ptr<HiggsGenEvent> gen( genEvt );
  genEvt->print();
  iEvent.put( gen );
}

// ------------ method called once each job just before starting event loop  ------------
void 
HiggsGenEventReco::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HiggsGenEventReco::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
HiggsGenEventReco::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
HiggsGenEventReco::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
HiggsGenEventReco::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
HiggsGenEventReco::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HiggsGenEventReco::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src");
  descriptions.add("genEvtHiggs", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HiggsGenEventReco);
