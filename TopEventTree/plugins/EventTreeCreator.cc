// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TopMass/TopEventTree/plugins/EventTreeCreator.h"

EventTreeCreator::EventTreeCreator(const edm::ParameterSet& cfg)
{
}

void
EventTreeCreator::analyze(const edm::Event& evt, const edm::EventSetup& setup)
{
}

void 
EventTreeCreator::beginJob()
{
  edm::Service<TFileService> fs;
  if( !fs ) throw edm::Exception( edm::errors::Configuration, "TFile Service is not registered in cfg file" );
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // book tree
  //////////////////////////////////////////////////////////////////////////////////////////////////
   
  eventTree = fs->make<TTree>("eventTree", "Tree for UHH top-quark analysis\nParticles are in order {TTBar, HadTop, LepTop, HadW, LepW, HadB, LightQ, LightQBar, LepB, Lepton, Neutrino}"); 
}

void
EventTreeCreator::endJob() 
{
}

EventTreeCreator::~EventTreeCreator()
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(EventTreeCreator);
