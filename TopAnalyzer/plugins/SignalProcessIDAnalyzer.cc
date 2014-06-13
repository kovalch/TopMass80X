#include "TopAnalysis/TopAnalyzer/plugins/SignalProcessIDAnalyzer.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

SignalProcessIDAnalyzer::SignalProcessIDAnalyzer(const ParameterSet& cfg)
{
}

SignalProcessIDAnalyzer::~SignalProcessIDAnalyzer()
{
}

void
SignalProcessIDAnalyzer::beginJob()
{
  edm::Service<TFileService> fs;
  if( !fs ){
    throw edm::Exception( edm::errors::Configuration,
                          "TFile Service is not registered in cfg file" );
  }

  id_= fs->make<TH1D>( "id", "id", 30, 0, 30);
  id_->GetXaxis()->SetTitle("SignalProcessID");
  id_->GetYaxis()->SetTitle("N");
}

void
SignalProcessIDAnalyzer::analyze(const Event& evt, const EventSetup&)
{
  edm::Handle<GenEventInfoProduct> GenEventInfoProduct;
  evt.getByLabel( "generator", GenEventInfoProduct );
  
  id_->Fill((*GenEventInfoProduct).signalProcessID());
}

void
SignalProcessIDAnalyzer::endJob()
{
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE( SignalProcessIDAnalyzer );
