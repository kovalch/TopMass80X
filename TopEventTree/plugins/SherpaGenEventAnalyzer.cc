/*
 * SherpaGenEventAnalyzer.cc
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
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TopMass/TopEventTree/plugins/SherpaGenEventAnalyzer.h"

SherpaGenEventAnalyzer::SherpaGenEventAnalyzer(const edm::ParameterSet& cfg):
  genJets_                (cfg.getParameter<edm::InputTag>("genJets"))
{
}

void
SherpaGenEventAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& setup)
{
  //////////////////////////////////////////////////////////////////////////////
  // INIT SherpaGenEvent
  ////////////////////////////////////////////////////////////////////////////

  sherpa->init();

  //////////////////////////////////////////////////////////////////////////
  // GENJETS
  ////////////////////////////////////////////////////////////////////////

  edm::Handle<std::vector< reco::GenJet > > genJets;
  evt.getByLabel(genJets_, genJets);

  //////////////////////////////////////////////////////////////////////////
  // PARTONS
  ////////////////////////////////////////////////////////////////////////
  
  edm::Handle<std::vector< reco::GenParticle > > genParticles;
  evt.getByLabel("genParticles", genParticles);
  
  for (std::vector< reco::GenParticle >::const_iterator p = genParticles->begin(); p != genParticles->end(); ++p) {
    // Select partons from top decay
    if (p->pt() == 0. || p->eta() > 5. || p->status() != 3) continue;
    
    // get pdgid, get 4-vector, match to genjet
    sherpa->flavour.push_back(p->pdgId());
    sherpa->parton.push_back(TLorentzVector(p->px(), p->py(), p->pz(), p->energy()));
    
    bool foundMatch = false;
    for (std::vector< reco::GenJet >::const_iterator ijet = genJets->begin(); ijet != genJets->end(); ++ijet) {
      double deta = p->eta() - ijet->eta();
      double dphi = TVector2::Phi_mpi_pi(p->phi() - ijet->phi());
      double dr   = sqrt( deta*deta + dphi*dphi );
      
      // Simple dR match of parton and GenJet
      // TODO: Save dR
      if (dr < 0.5) {
        sherpa->genJet.push_back(TLorentzVector(ijet->px(), ijet->py(), ijet->pz(), ijet->energy()));
        sherpa->dr.pushback(dr);
        foundMatch = true;
        break;
      }
    }
    if (!foundMatch) {
      sherpa->genJet.push_back(TLorentzVector(0,0,0,0));
      sherpa->dr.push_back(-1.);
    }
  }
  
  trs->Fill();
}

void
SherpaGenEventAnalyzer::beginJob()
{
  if( !trs ) throw edm::Exception( edm::errors::Configuration, "TreeRegistryService is not registered in cfg file!" );

  sherpa = new SherpaGenEvent();
  trs->Branch("sherpa.", sherpa);
}

void
SherpaGenEventAnalyzer::endJob()
{
}

SherpaGenEventAnalyzer::~SherpaGenEventAnalyzer()
{
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SherpaGenEventAnalyzer);
