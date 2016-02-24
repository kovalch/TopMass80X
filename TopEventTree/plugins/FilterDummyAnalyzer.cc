// -*- C++ -*-
//
// Package:    TopMass/FilterDummyAnalyzer
// Class:      FilterDummyAnalyzer
// 
/**\class FilterDummyAnalyzer FilterDummyAnalyzer.cc TopMass/FilterDummyAnalyzer/plugins/FilterDummyAnalyzer.cc

 Description: A Analyzer for the "FilterDummy"

 Implementation:
     Labels in "CutFlags.cc" and of the Filterprocesses have to be the same!
*/
//
// Original Author:  Christoph Garbers
//         Created:  Tue, 20 Oct 2015 09:51:54 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <string>
#include <vector>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TopMass/TopEventTree/interface/CutFlags.h"
#include "TopMass/TopEventTree/interface/TreeRegistryService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

//
// class declaration
//

class FilterDummyAnalyzer : public edm::EDAnalyzer {
   public:
      explicit FilterDummyAnalyzer(const edm::ParameterSet&);
      ~FilterDummyAnalyzer();

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      void setCutFlags(const edm::Event&);

      // ----------member data ---------------------------

	std::vector<edm::EDGetTokenT<bool> > vTokens_; //FIXME...

	edm::Service<TreeRegistryService> trs;
	CutFlags* flags_;
  	std::string cutFlagsBranchName_;

};

//
// constructors and destructor
//
FilterDummyAnalyzer::FilterDummyAnalyzer(const edm::ParameterSet& iConfig)
{
  flags_ = new CutFlags();
 for( unsigned int i = 0 ; i < flags_->getNumberOfFilterLabels() ; i++){
    vTokens_.push_back( mayConsume<bool>( "FilterDummy" + flags_->getFilterLabel(i) ) ); //FIXME
 }

   //now do what ever initialization is needed
   cutFlagsBranchName_ = iConfig.getParameter<std::string>("CutFlagsBranchName"); // set name of branch to be used to save
   flags_ = 0;
}


FilterDummyAnalyzer::~FilterDummyAnalyzer(){}

//
// member functions
//

// ------------ method called for each event  ------------
void
FilterDummyAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

  flags_->init();

  setCutFlags(iEvent);  

  trs->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
FilterDummyAnalyzer::beginJob()
{
   if( !trs ) throw edm::Exception( edm::errors::Configuration, "TreeRegistryService is not registered in cfg file!" );

  flags_ = new CutFlags();
  trs->Branch(cutFlagsBranchName_.c_str(), flags_);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
FilterDummyAnalyzer::endJob() {}


void
FilterDummyAnalyzer::setCutFlags(const edm::Event& iEvent)
{

using namespace edm;

 for( unsigned int i = 0 ; i < flags_->getNumberOfFilterLabels() ; i++){
	Handle<bool> hFilter;
	iEvent.getByLabel( "FilterDummy" + flags_->getFilterLabel(i) , hFilter );
 	if( !hFilter.failedToGet() && hFilter.isValid() )
		if( hFilter.product() )
			flags_->setFlag(1<<i);
 }
}

//define this as a plug-in
DEFINE_FWK_MODULE(FilterDummyAnalyzer);
