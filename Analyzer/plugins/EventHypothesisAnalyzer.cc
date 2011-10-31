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
#include <DataFormats/PatCandidates/interface/Muon.h>
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "LHAPDF/LHAPDF.h"

#include <TMath.h>

#include "TopMass/Analyzer/plugins/EventHypothesisAnalyzer.h"

EventHypothesisAnalyzer::EventHypothesisAnalyzer(const edm::ParameterSet& cfg):
  semiLepEvt_  (cfg.getParameter<edm::InputTag>("semiLepEvent")),
  hypoClassKey_(cfg.getParameter<edm::InputTag>("hypoClassKey")),
  jets_        (cfg.getParameter<edm::InputTag>("jets")),
	noPtEtaJets_ (cfg.getParameter<edm::InputTag>("noPtEtaJets")),
  leps_        (cfg.getParameter<edm::InputTag>("leps")),
  VertexSrc_   (cfg.getParameter<edm::InputTag>("VertexSrc")),
  PUWeightSrc_ (cfg.getParameter<edm::InputTag>("PUWeightSrc")),
	PUWeightUpSrc_  (cfg.getParameter<edm::InputTag>("PUWeightUpSrc")),
	PUWeightDownSrc_(cfg.getParameter<edm::InputTag>("PUWeightDownSrc")),
  bWeightSrc_  (cfg.getParameter<edm::InputTag>("bWeightSrc")),
  bWeightSrc_bTagSFUp_    (cfg.getParameter<edm::InputTag>("bWeightSrc_bTagSFUp")),
  bWeightSrc_bTagSFDown_  (cfg.getParameter<edm::InputTag>("bWeightSrc_bTagSFDown")),
  bWeightSrc_misTagSFUp_  (cfg.getParameter<edm::InputTag>("bWeightSrc_misTagSFUp")),
  bWeightSrc_misTagSFDown_(cfg.getParameter<edm::InputTag>("bWeightSrc_misTagSFDown")),
  muWeightSrc_ (cfg.getParameter<edm::InputTag>("muWeightSrc")),
  savePDFWeights_(cfg.getParameter<bool>("savePDFWeights" ))
{
}

void
EventHypothesisAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& setup)
{
  if(savePDFWeights_) {
    edm::Handle<GenEventInfoProduct> pdfstuff;
    evt.getByLabel("generator", pdfstuff);
    float Q   = pdfstuff->pdf()->scalePDF;
    int id1   = pdfstuff->pdf()->id.first;
    int id2   = pdfstuff->pdf()->id.second;
    double x1 = pdfstuff->pdf()->x.first;
    double x2 = pdfstuff->pdf()->x.second;
    double w0 = 0.;
    for(unsigned i=0; i <=44; ++i) {
      LHAPDF::usePDFMember(1,i);
      double xpdf1 = LHAPDF::xfx(1, x1, Q, id1);
      double xpdf2 = LHAPDF::xfx(1, x2, Q, id2);
      if(i<1)
	      w0 = xpdf1 * xpdf2;
      else
	      pdfWeights[i-1] = xpdf1 * xpdf2 / w0;
    }
  }

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
  jetMultiplicity = jets.isValid() ? jets->size() : -1;
	
	/*
	std::cout << "Selected jets" << std::endl;
	for (int i = 0; i < jetMultiplicity; i++) {
		std::cout << jets->at(i).pt() << "\t" << jets->at(i).eta() << "\t" << jets->at(i).phi() << std::endl; 
	}
	//*/
	
	edm::Handle<std::vector<pat::Jet> > noPtEtaJets;
  evt.getByLabel(noPtEtaJets_, noPtEtaJets);
	noPtEtaJetMultiplicity = noPtEtaJets.isValid() ? noPtEtaJets->size() : -1;
	
	/*
	std::cout << "All jets" << std::endl;
	for (int i = 0; i < noPtEtaJetMultiplicity; i++) {
		std::cout << noPtEtaJets->at(i).pt() << "\t" << noPtEtaJets->at(i).eta() << "\t" << noPtEtaJets->at(i).phi() << std::endl; 
	}
	//*/
	
	nlJetPt = 0;
	
  for (int i = 0; i < noPtEtaJetMultiplicity; i++) {
    int match = -1;
    for (int j = 0; j < 4; j++) {
      if (abs(noPtEtaJets->at(i).pt() - jets->at(j).pt()) < 0.01) {
        //std::cout << "match = " << j << std::endl;
        match = j;
      }			
    }
    if (match == -1) {
      //std::cout << "no match" << std::endl;
      nlJetPt  = noPtEtaJets->at(i).pt();
      nlJetEta = noPtEtaJets->at(i).eta();
      //std::cout << "Found NLJet " << i << std::endl;
      //std::cout << noPtEtaJets->at(i).pt() << "\t" << noPtEtaJets->at(i).eta() << "\t" << noPtEtaJets->at(i).phi() << std::endl;
      break;
    }
  }
	
	noPtEtaJetPt = 0;
	leadingJetPt = 0;
	
  for (int i = 0; i < noPtEtaJetMultiplicity; i++) {
    int match = -1;
    for (int j = 0; j < 4; j++) {
      if (abs(noPtEtaJets->at(i).pt() - jets->at(j).pt()) < 0.01) {
        match = j;
      }			
    }
    if (match == -1) {
      noPtEtaJetPt += noPtEtaJets->at(i).pt();
    }
    else leadingJetPt += noPtEtaJets->at(i).pt();
  }
  
  edm::Handle<std::vector<pat::Muon> > leps;
  evt.getByLabel(leps_, leps);
  
  edm::Handle<std::vector<reco::Vertex> > vertecies_h;
  evt.getByLabel(VertexSrc_, vertecies_h);
  nVertex = vertecies_h.isValid() ? vertecies_h->size() : -1;
  
  edm::Handle<double> PUWeightSrc_h;
  evt.getByLabel(PUWeightSrc_, PUWeightSrc_h);
  PUWeight = PUWeightSrc_h.isValid() ? *PUWeightSrc_h : -100.;
	
	edm::Handle<double> PUWeightUpSrc_h;
  evt.getByLabel(PUWeightUpSrc_, PUWeightUpSrc_h);
  PUWeightUp = PUWeightUpSrc_h.isValid() ? *PUWeightUpSrc_h : -100.;
	
	edm::Handle<double> PUWeightDownSrc_h;
  evt.getByLabel(PUWeightDownSrc_, PUWeightDownSrc_h);
  PUWeightDown = PUWeightDownSrc_h.isValid() ? *PUWeightDownSrc_h : -100.;
  
  edm::Handle<double> bWeightSrc_h;
  evt.getByLabel(bWeightSrc_, bWeightSrc_h);
  bWeight = bWeightSrc_h.isValid() ? *bWeightSrc_h : -100.;
  
  edm::Handle<double> bWeightSrc_bTagSFUp_h;
  evt.getByLabel(bWeightSrc_bTagSFUp_, bWeightSrc_bTagSFUp_h);
  bWeight_bTagSFUp = bWeightSrc_bTagSFUp_h.isValid() ? *bWeightSrc_bTagSFUp_h : -100.;
  
  edm::Handle<double> bWeightSrc_bTagSFDown_h;
  evt.getByLabel(bWeightSrc_bTagSFDown_, bWeightSrc_bTagSFDown_h);
  bWeight_bTagSFDown = bWeightSrc_bTagSFDown_h.isValid() ? *bWeightSrc_bTagSFDown_h : -100.;
  
  edm::Handle<double> bWeightSrc_misTagSFUp_h;
  evt.getByLabel(bWeightSrc_misTagSFUp_, bWeightSrc_misTagSFUp_h);
  bWeight_misTagSFUp = bWeightSrc_misTagSFUp_h.isValid() ? *bWeightSrc_misTagSFUp_h : -100.;
  
  edm::Handle<double> bWeightSrc_misTagSFDown_h;
  evt.getByLabel(bWeightSrc_misTagSFDown_, bWeightSrc_misTagSFDown_h);
  bWeight_misTagSFDown = bWeightSrc_misTagSFDown_h.isValid() ? *bWeightSrc_misTagSFDown_h : -100.;
  
  edm::Handle<double> muWeightSrc_h;
  evt.getByLabel(muWeightSrc_, muWeightSrc_h);
  muWeight = muWeightSrc_h.isValid() ? *muWeightSrc_h : -100.;
  
  MCWeight = PUWeightSrc_h.isValid() ? PUWeight * bWeight * muWeight : 1;

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
    
    const reco::Candidate* lepTop      = semiLepEvt->leptonicDecayTop       (hypoClassKey, h);
    const reco::Candidate* lepTopRaw   = semiLepEvt->leptonicDecayTop       (hypoClassKeyMVA, hMVA);
    const reco::Candidate* lepW        = semiLepEvt->leptonicDecayW         (hypoClassKey, h);
    const reco::Candidate* lepWRaw     = semiLepEvt->leptonicDecayW         (hypoClassKeyMVA, hMVA);
    const reco::Candidate* lepB        = semiLepEvt->leptonicDecayB         (hypoClassKey, h);
    const reco::Candidate* lepBRaw     = semiLepEvt->leptonicDecayB         (hypoClassKeyMVA, hMVA);
    const reco::Candidate* lepton      = semiLepEvt->singleLepton           (hypoClassKey, h);
    const reco::Candidate* leptonRaw   = semiLepEvt->singleLepton           (hypoClassKeyMVA, hMVA);
    const reco::Candidate* neutrino    = semiLepEvt->singleNeutrino         (hypoClassKey, h);
    const reco::Candidate* neutrinoRaw = semiLepEvt->singleNeutrino         (hypoClassKeyMVA, hMVA);

    
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // get genParticles
    //////////////////////////////////////////////////////////////////////////////////////////////////

    const reco::Candidate* genHadTop = semiLepEvt->hadronicDecayTop();
    const reco::Candidate* genHadW   = semiLepEvt->hadronicDecayW();
    const reco::Candidate* genHadQ   = semiLepEvt->hadronicDecayQuark();
    const reco::Candidate* genHadB   = semiLepEvt->hadronicDecayB();

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // fill eventTree with pt, eta and the masses of the reconstructed particles
    //////////////////////////////////////////////////////////////////////////////////////////////////
    
    hadQPt     = hadQ->pt();
    hadQEta    = hadQ->eta();
    hadQMass   = hadQ->mass();
    hadQE      = hadQ->energy();
    hadQBSSV   = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LightQ]).bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
    hadQJC     = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LightQ]).jetCharge();
    
    hadQRawPt  = hadQRaw->pt();
    
    if (genHadQ) {
      genHadQPt  = genHadQ->pt();
    }
    
    hadQBarPt     = hadQBar->pt();
    hadQBarEta    = hadQBar->eta();
    hadQBarMass   = hadQBar->mass();
    hadQBarE      = hadQBar->energy();
    hadQBarBSSV   = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LightQBar]).bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
    hadQBarJC     = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LightQBar]).jetCharge();
    
    hadQBarRawPt  = hadQBarRaw->pt();
    
    hadWPt     = hadW->pt();
    hadWEta    = hadW->eta();
    hadWMass   = hadW->mass();
    hadWE      = hadW->energy();
    
    hadWRawMass  = hadWRaw->mass();
    hadWRawPt    = hadWRaw->pt();
    
    /*
    double T     = ROOT::Math::VectorUtil::Angle(hadQRaw->polarP4(), hadQBarRaw->polarP4());
    double sigT  = sqrt(pow(jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LightQ]).resolTheta(), 2)
                      + pow(jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LightQBar]).resolTheta(), 2) );
    double E1    = hadQRawE;
    double sigE1 = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LightQ]).resolE();
    double E2    = hadQBarRawE;
    double sigE2 = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LightQBar]).resolE();
    
    hadWRawSigM  = sqrt(
                     pow(sin(T), 2) * E1 * E2 / (2 * (1-cos(T))) * pow(sigT, 2)
                   + (1-cos(T)) * E2 / (2 * E1) * pow(sigE1, 2)
                   + (1-cos(T)) * E1 / (2 * E2) * pow(sigE2, 2)
                   );
                   
                   //std::cout << hadWRawSigM << std::endl;
    */
  
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
    hadBJC     = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::HadB]).jetCharge();
    
    hadBRawPt  = hadBRaw->pt();
    
    leptonPt   = lepton->pt();
    leptonEta  = lepton->eta();
    leptonC    = leps->at(0).charge();
    
    leptonRawPt = leptonRaw->pt();
    
    neutrinoPt = neutrino->pt();
    neutrinoE  = neutrino->energy();
    
    neutrinoRawPt = neutrinoRaw->pt();
    
    lepWPt     = lepW->pt();
    lepWEta    = lepW->eta();
    lepWMass   = lepW->mass();
    lepWE      = lepW->energy();
    
    lepWRawMass  = lepWRaw->mass();
    
    lepBPt     = lepB->pt();
    lepBEta    = lepB->eta();
    lepBMass   = lepB->mass();
    lepBE      = lepB->energy();
    lepBBSSV   = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LepB]).bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
    lepBJC     = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LepB]).jetCharge();
    
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
    hadTopRawPt   = hadTopRaw->pt();
    hadTopRawEta  = hadTopRaw->eta();
  
    if (genHadTop) {
      genHadTopPt   = genHadTop->pt();
      genHadTopEta  = genHadTop->eta();
      genHadTopMass = genHadTop->mass();
    }
    
    lepTopPt   = lepTop->pt();
    lepTopEta  = lepTop->eta();
    lepTopMass = lepTop->mass();
    lepTopE    = lepTop->energy();
    
    lepTopRawMass = lepTopRaw->mass();
    
    deltaRHadQHadQBar     = ROOT::Math::VectorUtil::DeltaR(hadQ->polarP4(), hadQBar->polarP4());
    deltaThetaHadQHadQBar = ROOT::Math::VectorUtil::Angle(hadQ->polarP4(), hadQBar->polarP4());
    deltaRHadWHadB        = ROOT::Math::VectorUtil::DeltaR(hadW->polarP4(), hadB->polarP4());
    deltaThetaHadWHadB    = ROOT::Math::VectorUtil::Angle(hadW->polarP4(), hadB->polarP4());
    if (genHadW && genHadB) genDeltaThetaHadWHadB = ROOT::Math::VectorUtil::Angle(genHadW->polarP4(), genHadB->polarP4());
    deltaRLepBLepton      = ROOT::Math::VectorUtil::DeltaR(lepton->polarP4(), lepB->polarP4());
    deltaThetaLepBLepton  = ROOT::Math::VectorUtil::Angle(lepton->polarP4(), lepB->polarP4());
    deltaRHadBLepB        = ROOT::Math::VectorUtil::DeltaR(hadB->polarP4(), lepB->polarP4());;
    deltaThetaHadBLepB    = ROOT::Math::VectorUtil::Angle(hadB->polarP4(), lepB->polarP4());;
    
    sumB   = hadB->polarP4() + lepB->polarP4();
    sumBPt = sumB.pt();
    
    TTBar     = hadTop->polarP4() + lepTop->polarP4();
    TTBarPt   = TTBar.pt();
    TTBarMass = TTBar.mass();
    
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
    
    double c = fabs(leptonC + lepBJC);
    cProb = TMath::Max(5.35127e-04 + 9.94457e-02 * c - 1.42802e-01 * c*c + 7.77414e-02 * c*c*c - 1.54590e-02 * c*c*c*c, 0.);
    
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
  eventTree->Branch("hadQJC", &hadQJC, "hadQJC/D");
  
  eventTree->Branch("hadQRawPt", &hadQRawPt, "hadQRawPt/D");
  eventTree->Branch("genHadQPt", &genHadQPt, "genHadQPt/D");
  
  eventTree->Branch("hadQBarPt", &hadQBarPt, "hadQBarPt/D");
  eventTree->Branch("hadQBarEta", &hadQBarEta, "hadQBarEta/D");
  eventTree->Branch("hadQBarMass", &hadQBarMass, "hadQBarMass/D");
  eventTree->Branch("hadQBarE", &hadQBarE, "hadQBarE/D");
  eventTree->Branch("hadQBarBSSV", &hadQBarBSSV, "hadQBarBSSV/D");
  eventTree->Branch("hadQBarJC", &hadQBarJC, "hadQBarJC/D");
  
  eventTree->Branch("hadQBarRawPt", &hadQBarRawPt, "hadQBarRawPt/D");
  
  eventTree->Branch("hadWPt", &hadWPt, "hadWPt/D");
  eventTree->Branch("hadWEta", &hadWEta, "hadWEta/D");
  eventTree->Branch("hadWMass", &hadWMass, "hadWMass/D");
  eventTree->Branch("hadWE", &hadWE, "hadWE/D");
  
  eventTree->Branch("hadWRawMass", &hadWRawMass, "hadWRawMass/D");
  eventTree->Branch("hadWRawPt", &hadWRawPt, "hadWRawPt/D");
  
  eventTree->Branch("genHadWPt", &genHadWPt, "genHadWPt/D");
  eventTree->Branch("genHadWEta", &genHadWEta, "genHadWEta/D");
  eventTree->Branch("genHadWMass", &genHadWMass, "genHadWMass/D");
  eventTree->Branch("genHadWE", &genHadWE, "genHadWE/D");
  
  eventTree->Branch("hadBPt", &hadBPt, "hadBPt/D");
  eventTree->Branch("hadBEta", &hadBEta, "hadBEta/D");
  eventTree->Branch("hadBMass", &hadBMass, "hadBMass/D");
  eventTree->Branch("hadBE", &hadBE, "hadBE/D");
  eventTree->Branch("hadBBSSV", &hadBBSSV, "hadBBSSV/D");
  eventTree->Branch("hadBJC", &hadBJC, "hadBJC/D");
  
  eventTree->Branch("hadBRawPt", &hadBRawPt, "hadBRawPt/D");
  
  eventTree->Branch("leptonPt", &leptonPt, "leptonPt/D");
  eventTree->Branch("leptonEta", &leptonEta, "leptonEta/D");
  eventTree->Branch("leptonC", &leptonC, "leptonC/D");
  
  eventTree->Branch("leptonRawPt", &leptonRawPt, "leptonRawPt/D");
  
  eventTree->Branch("neutrinoPt", &neutrinoPt, "neutrinoPt/D");
  eventTree->Branch("neutrinoE", &neutrinoE, "neutrinoE/D");
  
  eventTree->Branch("neutrinoRawPt", &neutrinoRawPt, "neutrinoRawPt/D");
  
  eventTree->Branch("lepWPt", &lepWPt, "lepWPt/D");
  eventTree->Branch("lepWEta", &lepWEta, "lepWEta/D");
  eventTree->Branch("lepWMass", &lepWMass, "lepWMass/D");
  eventTree->Branch("lepWE", &lepWE, "lepWE/D");
  
  eventTree->Branch("lepWRawMass", &lepWRawMass, "lepWRawMass/D");
  
  eventTree->Branch("lepBPt", &lepBPt, "lepBPt/D");
  eventTree->Branch("lepBEta", &lepBEta, "lepBEta/D");
  eventTree->Branch("lepBMass", &lepBMass, "lepBMass/D");
  eventTree->Branch("lepBE", &lepBE, "lepBE/D");
  eventTree->Branch("lepBBSSV", &lepBBSSV, "lepBBSSV/D");
  eventTree->Branch("lepBJC", &lepBJC, "lepBJC/D");
  
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
  eventTree->Branch("hadTopRawPt", &hadTopRawPt, "hadTopRawPt/D");
  eventTree->Branch("hadTopRawEta", &hadTopRawEta, "hadTopRawEta/D");
  
  eventTree->Branch("lepTopPt", &lepTopPt, "lepTopPt/D");
  eventTree->Branch("lepTopEta", &lepTopEta, "lepTopEta/D");
  eventTree->Branch("lepTopMass", &lepTopMass, "lepTopMass/D");
  eventTree->Branch("lepTopE", &lepTopE, "lepTopE/D");
  
  eventTree->Branch("lepTopRawMass", &lepTopRawMass, "lepTopRawMass/D");
  
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
  eventTree->Branch("deltaRHadBLepB", &deltaRHadBLepB, "deltaRHadBLepB/D");
  eventTree->Branch("deltaThetaHadBLepB", &deltaThetaHadBLepB, "deltaThetaHadBLepB/D");
  
  eventTree->Branch("sumBPt", &sumBPt, "sumBPt/D");
  eventTree->Branch("TTBarPt", &TTBarPt, "TTBarPt/D");
  eventTree->Branch("TTBarMass", &TTBarMass, "TTBarMass/D");
  
  eventTree->Branch("jetMultiplicity", &jetMultiplicity, "jetMultiplicity/I");
	eventTree->Branch("noPtEtaJetMultiplicity", &noPtEtaJetMultiplicity, "noPtEtaJetMultiplicity/I");
	
	eventTree->Branch("noPtEtaJetPt", &noPtEtaJetPt, "noPtEtaJetPt/D");
	eventTree->Branch("leadingJetPt", &leadingJetPt, "leadingJetPt/D");
	
	eventTree->Branch("nlJetPt", &nlJetPt, "nlJetPt/D");
	eventTree->Branch("nlJetEta", &nlJetEta, "nlJetEta/D");
  
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
  eventTree->Branch("cProb", &cProb, "cProb/D");
  
  eventTree->Branch("nVertex", &nVertex, "nVertex/I");
  eventTree->Branch("PUWeight", &PUWeight, "PUWeight/D");
	eventTree->Branch("PUWeightUp", &PUWeightUp, "PUWeightUp/D");
	eventTree->Branch("PUWeightDown", &PUWeightDown, "PUWeightDown/D");
  eventTree->Branch("bWeight", &bWeight, "bWeight/D");
  eventTree->Branch("bWeight_bTagSFUp", &bWeight_bTagSFUp, "bWeight_bTagSFUp/D");
  eventTree->Branch("bWeight_bTagSFDown", &bWeight_bTagSFDown, "bWeight_bTagSFDown/D");
  eventTree->Branch("bWeight_misTagSFUp", &bWeight_misTagSFUp, "bWeight_misTagSFUp/D");
  eventTree->Branch("bWeight_misTagSFDown", &bWeight_misTagSFDown, "bWeight_misTagSFDown/D");
  eventTree->Branch("muWeight", &muWeight, "muWeight/D");
  eventTree->Branch("MCWeight", &MCWeight, "MCWeight/D");
  if(savePDFWeights_) {
    LHAPDF::initPDFSet(1, "cteq66.LHgrid");
    eventTree->Branch("pdfWeights", &pdfWeights, "pdfWeights[44]/D");
  }
  
  eventTree->Branch("target", &target, "target/I");

}

double EventHypothesisAnalyzer::QBTagProbabilitySSV(double bDiscriminator) {
  if (bDiscriminator < 1) return 0.7555;

  double p0 =  3.66166e-01;
  double p1 = -1.11745e+00;

  return exp(p0+p1*bDiscriminator);
}

double EventHypothesisAnalyzer::QBTagProbabilitySSVHEM(double bDiscriminator) {
  //TODO e(pt, eta) aus DB?
  double eb = 0.405;
  double el = 0.0084;
  
  if (bDiscriminator > 1.74) {
    return el/(eb+el);
  }
  else {
    return (1-el)/((1-eb)+(1-el));
  }
}

void
EventHypothesisAnalyzer::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(EventHypothesisAnalyzer);
