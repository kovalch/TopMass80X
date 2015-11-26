// -*- C++ -*-
//
// Package:    TopMass/FilterDummy
// Class:      FilterDummy
// 
/**\class FilterDummy FilterDummy.cc TopMass/FilterDummy/plugins/FilterDummy.cc

 Description: Dummy to count Filter passes

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Christoph Garbers
//         Created:  Mon, 19 Oct 2015 15:43:50 GMT
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
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//
// class declaration
//

class FilterDummy : public edm::EDProducer {
   public:
      explicit FilterDummy(const edm::ParameterSet&);
      ~FilterDummy();	

   private:
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;      
};

FilterDummy::FilterDummy(const edm::ParameterSet& iConfig)
{
   produces<bool>();  
}


FilterDummy::~FilterDummy()
{
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
FilterDummy::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

	std::auto_ptr<bool> fDoutput(new bool);
	*fDoutput=true;

	iEvent.put(fDoutput); 
}

// ------------ method called once each job just before starting event loop  ------------
void 
FilterDummy::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
FilterDummy::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(FilterDummy);
