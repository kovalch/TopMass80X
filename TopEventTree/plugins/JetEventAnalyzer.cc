/*
 * JetEventAnalyzer.cc
 *
 *  Created on: Feb 6, 2013
 *      Author: eschliec
 */

//#include <memory>

// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TopMass/TopEventTree/interface/TreeRegistryService.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

#include "TopMass/TopEventTree/plugins/JetEventAnalyzer.h"

JetEventAnalyzer::JetEventAnalyzer(const edm::ParameterSet& cfg):
jets_        (cfg.getParameter<edm::InputTag>("jets")),
//allJets_     (cfg.getParameter<edm::InputTag>("allJets")),
//noPtEtaJets_ (cfg.getParameter<edm::InputTag>("noPtEtaJets")),
kJetMAX_(cfg.getParameter<int>("maxNJets"))
{
}

void
JetEventAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& setup)
{
  //////////////////////////////////////////////////////////////////////////////
  // INIT JetEvent
  ////////////////////////////////////////////////////////////////////////////

  jet->init();

  //////////////////////////////////////////////////////////////////////////
  // JETS
  ////////////////////////////////////////////////////////////////////////

  edm::Handle<std::vector<pat::Jet> > jets;
  evt.getByLabel(jets_, jets);

  for(std::vector< pat::Jet >::const_iterator ijet = jets->begin(); ijet != jets->end(); ++ijet) {
    // write only kJetMAX_ jets into the event
    if(ijet - jets->begin() == kJetMAX_) break;

    jet->jet.push_back(TLorentzVector(ijet->px(), ijet->py(), ijet->pz(), ijet->energy()));

    jet->jetCharge.push_back(ijet->jetCharge());
    jet->jetFlavour.push_back(ijet->partonFlavour());
    jet->jetCSV.push_back(ijet->bDiscriminator("combinedSecondaryVertexBJetTags"));
  }

  trs->Fill();
}

void
JetEventAnalyzer::beginJob()
{
  if( !trs ) throw edm::Exception( edm::errors::Configuration, "TreeRegistryService is not registered in cfg file!" );

  jet = new JetEvent();
  trs->Branch("jet", jet);
}

void
JetEventAnalyzer::endJob()
{
}

JetEventAnalyzer::~JetEventAnalyzer()
{
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(JetEventAnalyzer);
