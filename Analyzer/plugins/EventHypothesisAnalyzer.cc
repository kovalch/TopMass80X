#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "AnalysisDataFormats/TopObjects/interface/TtSemiLeptonicEvent.h"
#include "AnalysisDataFormats/TopObjects/interface/TtSemiLepEvtPartons.h"
#include <DataFormats/PatCandidates/interface/Jet.h>

#include <TMath.h>

#include "TopMass/Analyzer/plugins/EventHypothesisAnalyzer.h"

EventHypothesisAnalyzer::EventHypothesisAnalyzer(const edm::ParameterSet& cfg):
  semiLepEvt_  (cfg.getParameter<edm::InputTag>("semiLepEvent")),
  hypoClassKey_(cfg.getParameter<edm::InputTag>("hypoClassKey"))
{
}

void
EventHypothesisAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& setup)
{

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // get a handle for the TtSemiLeptonicEvent and a key to the hypothesis
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  edm::Handle<TtSemiLeptonicEvent> semiLepEvt;
  evt.getByLabel(semiLepEvt_, semiLepEvt);

  edm::Handle<int> hypoClassKeyHandle;
  evt.getByLabel(hypoClassKey_, hypoClassKeyHandle);
  TtSemiLeptonicEvent::HypoClassKey& hypoClassKey = (TtSemiLeptonicEvent::HypoClassKey&) *hypoClassKeyHandle;
  
  edm::Handle<std::vector<pat::Jet> > jets;
  evt.getByLabel("selectedPatJetsAK5PF", jets);

  for(unsigned h=0; h<semiLepEvt->numberOfAvailableHypos(hypoClassKey); h++) {

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // check if hypothesis is available and valid in this event
    //////////////////////////////////////////////////////////////////////////////////////////////////

    if( !semiLepEvt->isHypoValid(hypoClassKey, h) ){
      edm::LogInfo("EventHypothesisAnalyzer") << "Hypothesis not valid for this event";
      return;
    }
    
    std::vector<int> jetLeptonCombinationCurrent = semiLepEvt->jetLeptonCombination(hypoClassKey, h);
    
    if (semiLepEvt->isHypoValid("kGenMatch") ) {
      std::vector<int> jetLeptonCombinationCurrent2 = semiLepEvt->jetLeptonCombination(hypoClassKey, h);
      std::swap(jetLeptonCombinationCurrent2[0], jetLeptonCombinationCurrent2[1]);
      std::vector<int> jetLeptonCombinationCurrentSorted = semiLepEvt->jetLeptonCombination(hypoClassKey, h);
      std::sort(jetLeptonCombinationCurrentSorted.begin(), jetLeptonCombinationCurrentSorted.begin()+4);

      std::vector<int> jetLeptonCombinationGenMatch = semiLepEvt->jetLeptonCombination("kGenMatch");
      std::vector<int> jetLeptonCombinationGenMatchSorted = semiLepEvt->jetLeptonCombination("kGenMatch");
      std::sort(jetLeptonCombinationGenMatchSorted.begin(), jetLeptonCombinationGenMatchSorted.begin()+4);

      std::vector<int> intersection(4);
      std::vector<int>::iterator intersection_it;

      intersection_it = std::set_intersection(jetLeptonCombinationCurrentSorted.begin(), jetLeptonCombinationCurrentSorted.begin()+4, jetLeptonCombinationGenMatchSorted.begin(), jetLeptonCombinationGenMatchSorted.begin()+4, intersection.begin());

      int maxNJets = 4;
      int maxMatchedJet = *max_element(jetLeptonCombinationGenMatch.begin(), jetLeptonCombinationGenMatch.begin()+4);

      // missing jet
      if (maxMatchedJet >= maxNJets) {
	target = -2;
      }
      // wrong jets
      else if (int(intersection_it - intersection.begin()) != 4) {
	target = -1;
      }
      // correct permutation
      else if (jetLeptonCombinationCurrent  == jetLeptonCombinationGenMatch ||
          jetLeptonCombinationCurrent2 == jetLeptonCombinationGenMatch) {
        target = 1;
      }
      // wrong permutation of correct jets
      else {
        target = 0;
      }
      genMatchDr = semiLepEvt->genMatchSumDR(h);
    }
    else {
      target = -10;
    }
    //std::cout << "\ttarget: " << target << " (" << hypoClassKey << ")" << std::endl;

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // get event Id
    //////////////////////////////////////////////////////////////////////////////////////////////////

    run             = evt.eventAuxiliary().run();
    luminosityBlock = evt.eventAuxiliary().luminosityBlock();
    event           = evt.eventAuxiliary().event();
    combi           = h;

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // get reconstructed particles from the hypothesis
    //////////////////////////////////////////////////////////////////////////////////////////////////

    const reco::Candidate* hadTop   = semiLepEvt->hadronicDecayTop       (hypoClassKey, h);
    const reco::Candidate* hadW     = semiLepEvt->hadronicDecayW         (hypoClassKey, h);
    const reco::Candidate* hadB     = semiLepEvt->hadronicDecayB         (hypoClassKey, h);
    const reco::Candidate* hadQ     = semiLepEvt->hadronicDecayQuark     (hypoClassKey, h);
    const reco::Candidate* hadQbar  = semiLepEvt->hadronicDecayQuarkBar  (hypoClassKey, h);
    
    //const reco::Candidate* lepTop   = semiLepEvt->leptonicDecayTop       (hypoClassKey, h);
    //const reco::Candidate* lepW     = semiLepEvt->leptonicDecayW         (hypoClassKey, h);
    const reco::Candidate* lepB     = semiLepEvt->leptonicDecayB         (hypoClassKey, h);
    const reco::Candidate* lepton   = semiLepEvt->singleLepton           (hypoClassKey, h);
    //const reco::Candidate* neutrino = semiLepEvt->singleNeutrino         (hypoClassKey, h);

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // get genParticles
    //////////////////////////////////////////////////////////////////////////////////////////////////

    const reco::Candidate* genHadTop = semiLepEvt->hadronicDecayTop();
    const reco::Candidate* genHadW   = semiLepEvt->hadronicDecayW();
    const reco::Candidate* genHadB   = semiLepEvt->hadronicDecayB();

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // fill eventTree with pt, eta and the masses of the reconstructed particles
    //////////////////////////////////////////////////////////////////////////////////////////////////
  
    hadQPt     = hadQ->pt();
    hadQEta    = hadQ->eta();
    hadQMass   = hadQ->mass();
    hadQE      = hadQ->energy();
    hadQBTCHE  = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LightQ]).bDiscriminator("trackCountingHighEffBJetTags");
    hadQBVMVA  = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LightQ]).bDiscriminator("combinedSecondaryVertexMVABJetTags");

    hadQbarPt     = hadQbar->pt();
    hadQbarEta    = hadQbar->eta();
    hadQbarMass   = hadQbar->mass();
    hadQbarE      = hadQbar->energy();
    hadQbarBTCHE  = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LightQBar]).bDiscriminator("trackCountingHighEffBJetTags");
    hadQbarBVMVA  = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LightQBar]).bDiscriminator("combinedSecondaryVertexMVABJetTags");
    
    hadWPt     = hadW->pt();
    hadWEta    = hadW->eta();
    hadWMass   = hadW->mass();
    hadWE      = hadW->energy();
  
    if (genHadW) {
      genHadWPt     = genHadW->pt();
      genHadWEta    = genHadW->eta();
      genHadWMass   = genHadW->mass();
      genHadWE      = genHadW->energy();
    }
  
    hadBPt     = hadB->pt();
    hadBEta    = hadB->eta();
    hadBMass   = hadB->mass();
    hadBE      = hadB->energy();
    hadBBTCHE  = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::HadB]).bDiscriminator("trackCountingHighEffBJetTags");
    hadBBVMVA  = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::HadB]).bDiscriminator("combinedSecondaryVertexMVABJetTags");
    
    lepBPt     = lepB->pt();
    lepBEta    = lepB->eta();
    lepBMass   = lepB->mass();
    lepBE      = lepB->energy();
    lepBBTCHE  = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LepB]).bDiscriminator("trackCountingHighEffBJetTags");
    lepBBVMVA  = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LepB]).bDiscriminator("combinedSecondaryVertexMVABJetTags");
  
    if (genHadB) {
      genHadBPt     = genHadB->pt();
      genHadBEta    = genHadB->eta();
      genHadBMass   = genHadB->mass();
      genHadBE      = genHadB->energy();
    }
  
    hadTopPt   = hadTop->pt();
    hadTopEta  = hadTop->eta();
    hadTopMass = hadTop->mass();
    hadTopE    = hadTop->energy();
    
    //std::cout << hadB->resolM() << std::endl;
    //std::cout << hadTop->resolM() << std::endl;
  
    if (genHadTop) {
      genHadTopPt   = genHadTop->pt();
      genHadTopEta  = genHadTop->eta();
      genHadTopMass = genHadTop->mass();
    }
    
    deltaRHadQHadQBar     = ROOT::Math::VectorUtil::DeltaR(hadQ->polarP4(), hadQbar->polarP4());
    deltaThetaHadQHadQBar = ROOT::Math::VectorUtil::Angle(hadQ->polarP4(), hadQbar->polarP4());
    deltaRHadWHadB        = ROOT::Math::VectorUtil::DeltaR(hadW->polarP4(), hadB->polarP4());
    deltaThetaHadWHadB    = ROOT::Math::VectorUtil::Angle(hadW->polarP4(), hadB->polarP4());
    if (genHadW && genHadB) genDeltaThetaHadWHadB = ROOT::Math::VectorUtil::Angle(genHadW->polarP4(), genHadB->polarP4());
    deltaRLepBLepton      = ROOT::Math::VectorUtil::DeltaR(lepton->polarP4(), lepB->polarP4());
    deltaThetaLepBLepton  = ROOT::Math::VectorUtil::Angle(lepton->polarP4(), lepB->polarP4());
    
    mvaDisc    = semiLepEvt->mvaDisc(h);
    fitChi2    = semiLepEvt->fitChi2(h);
    fitProb    = semiLepEvt->fitProb(h);
    hitFitChi2 = semiLepEvt->hitFitChi2(h);
    hitFitProb = semiLepEvt->hitFitProb(h);
    hitFitMT   = semiLepEvt->hitFitMT(h);
    hitFitSigMT= semiLepEvt->hitFitSigMT(h);
    bProb      = QBTagProbability(hadQBTCHE) * QBTagProbability(hadQbarBTCHE)
                 * (1 - QBTagProbability(hadBBTCHE))
                 * (1 - QBTagProbability(lepBBTCHE));
    hadBProb   = QBTagProbability(hadQBTCHE) * QBTagProbability(hadQbarBTCHE)
                 * (1 - QBTagProbability(hadBBTCHE));
  
    eventTree -> Fill();
  
  }

}

void 
EventHypothesisAnalyzer::beginJob()
{
  edm::Service<TFileService> fs;
  if( !fs ) throw edm::Exception( edm::errors::Configuration, "TFile Service is not registered in cfg file" );

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // book tree
  //////////////////////////////////////////////////////////////////////////////////////////////////
					       
  eventTree = fs->make<TTree>("eventTree", "Tree for kinematic variables"); 

  eventTree->Branch("run", &run, "run/I");
  eventTree->Branch("luminosityBlock", &luminosityBlock, "luminosityBlock/I");
  eventTree->Branch("event", &event, "event/I");
  eventTree->Branch("combi", &combi, "combi/I");
      
  eventTree->Branch("hadQPt", &hadQPt, "hadQPt/D");
  eventTree->Branch("hadQEta", &hadQEta, "hadQEta/D");
  eventTree->Branch("hadQMass", &hadQMass, "hadQMass/D");
  eventTree->Branch("hadQE", &hadQE, "hadQE/D");
  eventTree->Branch("hadQBTCHE", &hadQBTCHE, "hadQBTCHE/D");
  eventTree->Branch("hadQBVMVA", &hadQBVMVA, "hadQBVMVA/D");
  
  eventTree->Branch("hadQbarPt", &hadQbarPt, "hadQbarPt/D");
  eventTree->Branch("hadQbarEta", &hadQbarEta, "hadQbarEta/D");
  eventTree->Branch("hadQbarMass", &hadQbarMass, "hadQbarMass/D");
  eventTree->Branch("hadQbarE", &hadQbarE, "hadQbarE/D");
  eventTree->Branch("hadQbarBTCHE", &hadQbarBTCHE, "hadQbarBTCHE/D");
  eventTree->Branch("hadQbarBVMVA", &hadQbarBVMVA, "hadQbarBVMVA/D");
  
  eventTree->Branch("hadWPt", &hadWPt, "hadWPt/D");
  eventTree->Branch("hadWEta", &hadWEta, "hadWEta/D");
  eventTree->Branch("hadWMass", &hadWMass, "hadWMass/D");
  eventTree->Branch("hadWE", &hadWE, "hadWE/D");
  
  eventTree->Branch("genHadWPt", &genHadWPt, "genHadWPt/D");
  eventTree->Branch("genHadWEta", &genHadWEta, "genHadWEta/D");
  eventTree->Branch("genHadWMass", &genHadWMass, "genHadWMass/D");
  eventTree->Branch("genHadWE", &genHadWE, "genHadWE/D");
  
  eventTree->Branch("hadBPt", &hadBPt, "hadBPt/D");
  eventTree->Branch("hadBEta", &hadBEta, "hadBEta/D");
  eventTree->Branch("hadBMass", &hadBMass, "hadBMass/D");
  eventTree->Branch("hadBE", &hadBE, "hadBE/D");
  eventTree->Branch("hadBBTCHE", &hadBBTCHE, "hadBBTCHE/D");
  eventTree->Branch("hadBBVMVA", &hadBBVMVA, "hadBBVMVA/D");
  
  eventTree->Branch("lepBPt", &lepBPt, "lepBPt/D");
  eventTree->Branch("lepBEta", &lepBEta, "lepBEta/D");
  eventTree->Branch("lepBMass", &lepBMass, "lepBMass/D");
  eventTree->Branch("lepBE", &lepBE, "lepBE/D");
  eventTree->Branch("lepBBTCHE", &lepBBTCHE, "lepBBTCHE/D");
  eventTree->Branch("lepBBVMVA", &lepBBVMVA, "lepBBVMVA/D");
  
  eventTree->Branch("genHadBPt", &genHadBPt, "genHadBPt/D");
  eventTree->Branch("genHadBEta", &genHadBEta, "genHadBEta/D");
  eventTree->Branch("genHadBMass", &genHadBMass, "genHadBMass/D");
  eventTree->Branch("genHadBE", &genHadBE, "genHadBE/D");
  
  eventTree->Branch("hadTopPt", &hadTopPt, "hadTopPt/D");
  eventTree->Branch("hadTopEta", &hadTopEta, "hadTopEta/D");
  eventTree->Branch("hadTopMass", &hadTopMass, "hadTopMass/D");
  eventTree->Branch("hadTopE", &hadTopE, "hadTopE/D");
  
  eventTree->Branch("genHadTopPt", &genHadTopPt, "genHadTopPt/D");
  eventTree->Branch("genHadTopEta", &genHadTopEta, "genHadTopEta/D");
  eventTree->Branch("genHadTopMass", &genHadTopMass, "genHadTopMass/D");
  
  eventTree->Branch("deltaRHadQHadQBar", &deltaRHadQHadQBar, "deltaRHadQHadQBar/D");
  eventTree->Branch("deltaThetaHadQHadQBar", &deltaThetaHadQHadQBar, "deltaThetaHadQHadQBar/D");
  eventTree->Branch("deltaRHadWHadB", &deltaRHadWHadB, "deltaRHadWHadB/D");
  eventTree->Branch("deltaThetaHadWHadB", &deltaThetaHadWHadB, "deltaThetaHadWHadB/D");
  eventTree->Branch("genDeltaThetaHadWHadB", &genDeltaThetaHadWHadB, "genDeltaThetaHadWHadB/D");
  eventTree->Branch("deltaRLepBLepton", &deltaRLepBLepton, "deltaRLepBLepton/D");
  eventTree->Branch("deltaThetaLepBLepton", &deltaThetaLepBLepton, "deltaThetaLepBLepton/D");
  
  eventTree->Branch("genMatchDr", &genMatchDr, "genMatchDr/D");
  eventTree->Branch("mvaDisc", &mvaDisc, "mvaDisc/D");
  eventTree->Branch("fitChi2", &fitChi2, "fitChi2/D");
  eventTree->Branch("fitProb", &fitProb, "fitProb/D");
  eventTree->Branch("hitFitChi2", &hitFitChi2, "hitFitChi2/D");
  eventTree->Branch("hitFitProb", &hitFitProb, "hitFitProb/D");
  eventTree->Branch("hitFitMT", &hitFitMT, "hitFitMT/D");
  eventTree->Branch("hitFitSigMT", &hitFitSigMT, "hitFitSigMT/D");
  eventTree->Branch("bProb", &bProb, "bProb/D");
  eventTree->Branch("hadBProb", &hadBProb, "hadBProb/D");
  
  eventTree->Branch("target", &target, "target/I");

}

double EventHypothesisAnalyzer::QBTagProbability(double bDiscriminator) {
  if (bDiscriminator == -100) return 0.787115;
  if (bDiscriminator < 0) return 1;
  
  double p0 = 5.91566e+00;
  double p1 = 5.94611e-01;
  double p2 = 3.53592e+00;
  
  return p0 * TMath::Voigt(bDiscriminator, p1, p2);
}

void
EventHypothesisAnalyzer::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(EventHypothesisAnalyzer);
