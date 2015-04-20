// -*- C++ -*-
//
// Package:    RunRangeFilter
// Class:      RunRangeFilter
//
/**\class RunRangeFilter RunRangeFilter.cc TopAnalysis/RunRangeFilter/src/RunRangeFilter.cc

Description: Filter for run-ranges to allow different trigger path for different runs.

Implementation:
Both limits are included in the valid range. If a limit is 0 it is not tested.*/
//
// Original Author:  Benjamin Lutz,,,DESY
//         Created:  Thu Dec 16 16:13:22 CET 2010
// $Id: RunRangeFilter.cc,v 1.1 2011/02/28 12:16:50 blutz Exp $
//
//


// system include files
#include <memory>
#include <iostream>

// user include files
//#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class RunRangeFilter : public edm::EDFilter {
public:
  explicit RunRangeFilter(const edm::ParameterSet&);
  ~RunRangeFilter();

private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  edm::RunNumber_t _minRunNumber;
  edm::RunNumber_t _maxRunNumber;

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
RunRangeFilter::RunRangeFilter(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed

  _minRunNumber = iConfig.getUntrackedParameter<edm::RunNumber_t>("minRunNumber",0);
  _maxRunNumber = iConfig.getUntrackedParameter<edm::RunNumber_t>("maxRunNumber",0);

}


RunRangeFilter::~RunRangeFilter()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
RunRangeFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  std::cout << "This filter is currently blid --> PLEASE FIX THIS TO WORK WITH CMSSW_3_11_X <--" << std::endl;
  /*
    PLEASE FIX THIS TO WORK WITH CMSSW_3_11_X 

    RunNumber_t runNumber = iEvent.getRun().id().run();
    
    if (_minRunNumber > 0 && runNumber < _minRunNumber ) return false;
    if (_maxRunNumber > 0 && runNumber > _maxRunNumber ) return false;
  */
  return true;
}

// ------------ method called once each job just before starting event loop  ------------
void
RunRangeFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
RunRangeFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(RunRangeFilter);
