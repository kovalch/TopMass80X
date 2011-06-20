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
  hypoClassKey_(cfg.getParameter<edm::InputTag>("hypoClassKey")),
  jets_        (cfg.getParameter<edm::InputTag>("jets"))
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
  
  TtSemiLeptonicEvent::HypoClassKey hypoClassKeyGenMatch = (TtSemiLeptonicEvent::HypoClassKey) 3;
  TtSemiLeptonicEvent::HypoClassKey hypoClassKeyMVA      = (TtSemiLeptonicEvent::HypoClassKey) 4;
  TtSemiLeptonicEvent::HypoClassKey hypoClassKeyKinFit   = (TtSemiLeptonicEvent::HypoClassKey) 5;
  TtSemiLeptonicEvent::HypoClassKey hypoClassKeyHitFit   = (TtSemiLeptonicEvent::HypoClassKey) 8;
  
  edm::Handle<std::vector<pat::Jet> > jets;
  evt.getByLabel(jets_, jets);

  for(unsigned h=0; h<semiLepEvt->numberOfAvailableHypos(hypoClassKey); h++) {

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // check if hypothesis is available and valid in this event
    //////////////////////////////////////////////////////////////////////////////////////////////////

    if( !semiLepEvt->isHypoValid(hypoClassKey, h) ){
      edm::LogInfo("EventHypothesisAnalyzer") << "Hypothesis not valid for this event";
      return;
    }
    
    int hGenMatch = semiLepEvt->correspondingHypo(hypoClassKey, h, hypoClassKeyGenMatch);
    int hMVA      = semiLepEvt->correspondingHypo(hypoClassKey, h, hypoClassKeyMVA);
    int hKinFit   = semiLepEvt->correspondingHypo(hypoClassKey, h, hypoClassKeyKinFit);
    int hHitFit;
    if (hypoClassKey == hypoClassKeyHitFit) hHitFit = h;
    else hHitFit  = semiLepEvt->correspondingHypo(hypoClassKey, h, hypoClassKeyHitFit);
    
    if( !semiLepEvt->isHypoValid(hypoClassKeyMVA, hMVA) ){
      edm::LogInfo("EventHypothesisAnalyzer") << "MVA Hypothesis not valid for this event";
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
      genMatchDr = semiLepEvt->genMatchSumDR(hGenMatch);
    }
    else {
      target = -10;
      genMatchDr = -10;
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
    
    const reco::Candidate* hadTop     = semiLepEvt->hadronicDecayTop       (hypoClassKey, h);
    const reco::Candidate* hadTopRaw  = semiLepEvt->hadronicDecayTop       (hypoClassKeyMVA, hMVA);
    const reco::Candidate* hadW       = semiLepEvt->hadronicDecayW         (hypoClassKey, h);
    const reco::Candidate* hadWRaw    = semiLepEvt->hadronicDecayW         (hypoClassKeyMVA, hMVA);
    const reco::Candidate* hadB       = semiLepEvt->hadronicDecayB         (hypoClassKey, h);
    const reco::Candidate* hadBRaw    = semiLepEvt->hadronicDecayB         (hypoClassKeyMVA, hMVA);
    const reco::Candidate* hadQ       = semiLepEvt->hadronicDecayQuark     (hypoClassKey, h);
    const reco::Candidate* hadQRaw    = semiLepEvt->hadronicDecayQuark     (hypoClassKeyMVA, hMVA);
    const reco::Candidate* hadQBar    = semiLepEvt->hadronicDecayQuarkBar  (hypoClassKey, h);
    const reco::Candidate* hadQBarRaw = semiLepEvt->hadronicDecayQuarkBar  (hypoClassKeyMVA, hMVA);
    
    //const reco::Candidate* lepTop   = semiLepEvt->leptonicDecayTop       (hypoClassKey, h);
    //const reco::Candidate* lepW     = semiLepEvt->leptonicDecayW         (hypoClassKey, h);
    const reco::Candidate* lepB     = semiLepEvt->leptonicDecayB         (hypoClassKey, h);
    const reco::Candidate* lepBRaw  = semiLepEvt->leptonicDecayB         (hypoClassKeyMVA, hMVA);
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
    hadQBSSV   = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LightQ]).bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
    
    hadQRawE   = hadQRaw->energy();
    
    hadQBarPt     = hadQBar->pt();
    hadQBarEta    = hadQBar->eta();
    hadQBarMass   = hadQBar->mass();
    hadQBarE      = hadQBar->energy();
    hadQBarBSSV   = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LightQBar]).bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
    
    hadQBarRawE   = hadQBarRaw->energy();
    
    hadWPt     = hadW->pt();
    hadWEta    = hadW->eta();
    hadWMass   = hadW->mass();
    hadWE      = hadW->energy();
    
    hadWRawMass   = hadWRaw->mass();
  
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
    hadBBSSV   = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::HadB]).bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
    
    hadBRawE   = hadBRaw->energy();
    
    lepBPt     = lepB->pt();
    lepBEta    = lepB->eta();
    lepBMass   = lepB->mass();
    lepBE      = lepB->energy();
    lepBBSSV   = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LepB]).bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
    
    lepBRawE   = lepBRaw->energy();
  
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
    
    hadTopRawMass = hadTopRaw->mass();
  
    if (genHadTop) {
      genHadTopPt   = genHadTop->pt();
      genHadTopEta  = genHadTop->eta();
      genHadTopMass = genHadTop->mass();
    }
    
    deltaRHadQHadQBar     = ROOT::Math::VectorUtil::DeltaR(hadQ->polarP4(), hadQBar->polarP4());
    deltaThetaHadQHadQBar = ROOT::Math::VectorUtil::Angle(hadQ->polarP4(), hadQBar->polarP4());
    deltaRHadWHadB        = ROOT::Math::VectorUtil::DeltaR(hadW->polarP4(), hadB->polarP4());
    deltaThetaHadWHadB    = ROOT::Math::VectorUtil::Angle(hadW->polarP4(), hadB->polarP4());
    if (genHadW && genHadB) genDeltaThetaHadWHadB = ROOT::Math::VectorUtil::Angle(genHadW->polarP4(), genHadB->polarP4());
    deltaRLepBLepton      = ROOT::Math::VectorUtil::DeltaR(lepton->polarP4(), lepB->polarP4());
    deltaThetaLepBLepton  = ROOT::Math::VectorUtil::Angle(lepton->polarP4(), lepB->polarP4());
    
    mvaDisc    = semiLepEvt->mvaDisc(hMVA);
    fitChi2    = semiLepEvt->fitChi2(hKinFit);
    fitProb    = semiLepEvt->fitProb(hKinFit);
    hitFitChi2 = semiLepEvt->hitFitChi2(hHitFit);
    hitFitProb = semiLepEvt->hitFitProb(hHitFit);
    hitFitMT   = semiLepEvt->hitFitMT(hHitFit);
    hitFitSigMT= semiLepEvt->hitFitSigMT(hHitFit);

    bProbSSV   = QBTagProbabilitySSV(hadQBSSV) * QBTagProbabilitySSV(hadQBarBSSV)
                 * (1 - QBTagProbabilitySSV(hadBBSSV))
                 * (1 - QBTagProbabilitySSV(lepBBSSV));
    hadBProbSSV= QBTagProbabilitySSV(hadQBSSV) * QBTagProbabilitySSV(hadQBarBSSV)
                 * (1 - QBTagProbabilitySSV(hadBBSSV));
  
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
  eventTree->Branch("hadQBSSV", &hadQBSSV, "hadQBSSV/D");
  
  eventTree->Branch("hadQRawE", &hadQRawE, "hadQRawE/D");
  
  eventTree->Branch("hadQBarPt", &hadQBarPt, "hadQBarPt/D");
  eventTree->Branch("hadQBarEta", &hadQBarEta, "hadQBarEta/D");
  eventTree->Branch("hadQBarMass", &hadQBarMass, "hadQBarMass/D");
  eventTree->Branch("hadQBarE", &hadQBarE, "hadQBarE/D");
  eventTree->Branch("hadQBarBSSV", &hadQBarBSSV, "hadQBarBSSV/D");
  
  eventTree->Branch("hadQBarRawE", &hadQBarRawE, "hadQBarRawE/D");
  
  eventTree->Branch("hadWPt", &hadWPt, "hadWPt/D");
  eventTree->Branch("hadWEta", &hadWEta, "hadWEta/D");
  eventTree->Branch("hadWMass", &hadWMass, "hadWMass/D");
  eventTree->Branch("hadWE", &hadWE, "hadWE/D");
  
  eventTree->Branch("hadWRawMass", &hadWRawMass, "hadWRawMass/D");
  
  eventTree->Branch("genHadWPt", &genHadWPt, "genHadWPt/D");
  eventTree->Branch("genHadWEta", &genHadWEta, "genHadWEta/D");
  eventTree->Branch("genHadWMass", &genHadWMass, "genHadWMass/D");
  eventTree->Branch("genHadWE", &genHadWE, "genHadWE/D");
  
  eventTree->Branch("hadBPt", &hadBPt, "hadBPt/D");
  eventTree->Branch("hadBEta", &hadBEta, "hadBEta/D");
  eventTree->Branch("hadBMass", &hadBMass, "hadBMass/D");
  eventTree->Branch("hadBE", &hadBE, "hadBE/D");
  eventTree->Branch("hadBBSSV", &hadBBSSV, "hadBBSSV/D");
  
  eventTree->Branch("hadBRawE", &hadBRawE, "hadBRawE/D");
  
  eventTree->Branch("lepBPt", &lepBPt, "lepBPt/D");
  eventTree->Branch("lepBEta", &lepBEta, "lepBEta/D");
  eventTree->Branch("lepBMass", &lepBMass, "lepBMass/D");
  eventTree->Branch("lepBE", &lepBE, "lepBE/D");
  eventTree->Branch("lepBBSSV", &lepBBSSV, "lepBBSSV/D");
  
  eventTree->Branch("lepBRawE", &lepBRawE, "lepBRawE/D");
  
  eventTree->Branch("genHadBPt", &genHadBPt, "genHadBPt/D");
  eventTree->Branch("genHadBEta", &genHadBEta, "genHadBEta/D");
  eventTree->Branch("genHadBMass", &genHadBMass, "genHadBMass/D");
  eventTree->Branch("genHadBE", &genHadBE, "genHadBE/D");
  
  eventTree->Branch("hadTopPt", &hadTopPt, "hadTopPt/D");
  eventTree->Branch("hadTopEta", &hadTopEta, "hadTopEta/D");
  eventTree->Branch("hadTopMass", &hadTopMass, "hadTopMass/D");
  eventTree->Branch("hadTopE", &hadTopE, "hadTopE/D");
  
  eventTree->Branch("hadTopRawMass", &hadTopRawMass, "hadTopRawMass/D");
  
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
  eventTree->Branch("bProbSSV", &bProbSSV, "bProbSSV/D");
  eventTree->Branch("hadBProbSSV", &hadBProbSSV, "hadBProbSSV/D");
  
  eventTree->Branch("target", &target, "target/I");

}

double EventHypothesisAnalyzer::QBTagProbabilitySSV(double bDiscriminator) {
  if (bDiscriminator < 1) return 0.7555;

  double p0 =  3.66166e-01;
  double p1 = -1.11745e+00;

  return exp(p0+p1*bDiscriminator);
}

void
EventHypothesisAnalyzer::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(EventHypothesisAnalyzer);
