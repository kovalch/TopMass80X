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
  
  // isr, b, bbar, wpq1, wpq2, wmq1, wmq2
  int partons[] = {-1, -1, -1, -1, -1, -1, -1};
  
  // Search for partons
  for (std::vector< reco::GenParticle >::const_iterator p = genParticles->begin(); p != genParticles->end(); ++p) {
    if ((p->pt() == 0.) || (fabs(p->eta()) > 5.)) continue;
    if (p->status() == 3) {
      if (abs(p->pdgId()) !=  6) partons[0] = p - genParticles->begin();
    }
    if ((p->status() == 11) && (abs(p->pdgId()) < 20 )) { // Only fermion daughters
      if ((p->mother()->pdgId() ==  6) && (p->pdgId() ==  5)) partons[1] = p - genParticles->begin();
      if ((p->mother()->pdgId() == -6) && (p->pdgId() == -5)) partons[2] = p - genParticles->begin();
      
      if ((p->mother()->mother()->pdgId() ==  6) && (p->mother()->pdgId() ==  24)) {
        if (p->pdgId() > 0) partons[3] = p - genParticles->begin();
        if (p->pdgId() < 0) partons[4] = p - genParticles->begin();
      }
      if ((p->mother()->mother()->pdgId() == -6) && (p->mother()->pdgId() == -24)) {
        if (p->pdgId() > 0) partons[5] = p - genParticles->begin();
        if (p->pdgId() < 0) partons[6] = p - genParticles->begin();
      }
    }
  }
  
  // Init top and W
  sherpa->genJetTop.push_back(TLorentzVector(0,0,0,0)); sherpa->genJetTop.push_back(TLorentzVector(0,0,0,0));
  sherpa->genJetW  .push_back(TLorentzVector(0,0,0,0)); sherpa->genJetW  .push_back(TLorentzVector(0,0,0,0));
  
  // Save partons and genJets
  for (int i = 0; i < 7; ++i) {
    if (partons[i] == -1) continue;
    reco::GenParticle p = genParticles->at(partons[i]);
    sherpa->flavour.push_back(p.pdgId());
    sherpa->parton.push_back(TLorentzVector(p.px(), p.py(), p.pz(), p.energy()));
    
    bool foundMatch = false;
    if (abs(p.pdgId()) < 6 || p.pdgId() == 21) {
      for (std::vector< reco::GenJet >::const_iterator ijet = genJets->begin(); ijet != genJets->end(); ++ijet) {
        double deta = p.eta() - ijet->eta();
        double dphi = TVector2::Phi_mpi_pi(p.phi() - ijet->phi());
        double dr   = sqrt( deta*deta + dphi*dphi );
        
        // Simple dR match of parton and GenJet
        if (dr < 0.5) {
          sherpa->genJet.push_back(TLorentzVector(ijet->px(), ijet->py(), ijet->pz(), ijet->energy()));
          sherpa->dr.push_back(dr);
          foundMatch = true;
          
          if (i==1 || i==3 || i==4) sherpa->genJetTop[0] += sherpa->genJet.back();
          if (        i==3 || i==4) sherpa->genJetW  [0] += sherpa->genJet.back();
          if (i==2 || i==5 || i==6) sherpa->genJetTop[1] += sherpa->genJet.back();
          if (        i==5 || i==6) sherpa->genJetW  [1] += sherpa->genJet.back();
          
          break;
        }
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
