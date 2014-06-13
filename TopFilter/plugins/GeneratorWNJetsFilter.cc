#include "TopAnalysis/TopFilter/plugins/GeneratorWNJetsFilter.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/Framework/interface/Event.h"

GeneratorWNJetsFilter::GeneratorWNJetsFilter(const edm::ParameterSet& cfg):
  NJet_( cfg.getParameter<int>("NJet") )
{
}

bool GeneratorWNJetsFilter::filter(edm::Event& event, const edm::EventSetup& setup)
{
  bool pass = false;

  edm::Handle<GenEventInfoProduct> GenEventInfoProduct;
  event.getByLabel( "generator", GenEventInfoProduct );

  if (((*GenEventInfoProduct).signalProcessID() == 0 && NJet_ == 0)
   || ((*GenEventInfoProduct).signalProcessID() == 1 && NJet_ == 1)
   || ((*GenEventInfoProduct).signalProcessID() == 2 && NJet_ == 2)
   || ((*GenEventInfoProduct).signalProcessID() == 3 && NJet_ == 3)
   || ((*GenEventInfoProduct).signalProcessID() == 4 && NJet_ == 4)
   || ((*GenEventInfoProduct).signalProcessID() == 5 && NJet_ == 0)
   || ((*GenEventInfoProduct).signalProcessID() == 6 && NJet_ == 1)
   || ((*GenEventInfoProduct).signalProcessID() == 7 && NJet_ == 2)
   || ((*GenEventInfoProduct).signalProcessID() == 8 && NJet_ == 3)
   || ((*GenEventInfoProduct).signalProcessID() == 9 && NJet_ == 4)) {
	
	  pass = true;
  }

  return pass;
}

void GeneratorWNJetsFilter::beginJob()
{
}

void GeneratorWNJetsFilter::endJob()
{
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE( GeneratorWNJetsFilter );
