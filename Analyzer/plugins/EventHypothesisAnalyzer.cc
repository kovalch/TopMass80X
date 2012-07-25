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
  allJets_     (cfg.getParameter<edm::InputTag>("allJets")),
  noPtEtaJets_ (cfg.getParameter<edm::InputTag>("noPtEtaJets")),
  leps_        (cfg.getParameter<edm::InputTag>("leps")),
  mets_        (cfg.getParameter<edm::InputTag>("mets")),

  PUSrc_       (cfg.getParameter<edm::InputTag>("PUSrc")),
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
  mcWeightSrc_ (cfg.getParameter<edm::InputTag>("mcWeightSrc")),
  
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
  
  for (int i = 0; i < std::min(jetMultiplicity, 6); ++i) {
    jetsPt[i]   = jets->at(i).pt();
    jetsEta[i]  = jets->at(i).eta();
    jetsPhi[i]  = jets->at(i).phi();
    jetsBCSV[i] = jets->at(i).bDiscriminator("combinedSecondaryVertexBJetTags");
  }
  
  edm::Handle<std::vector<pat::Jet> > allJets;
  evt.getByLabel(allJets_, allJets);
  allJetMultiplicity = allJets.isValid() ? allJets->size() : -1;

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

  edm::Handle<std::vector<pat::Jet> > bottomSSVJets;
  evt.getByLabel("tightBottomSSVPFJets", bottomSSVJets);
  bottomSSVJetMultiplicity = bottomSSVJets.isValid() ? bottomSSVJets->size() : -1;

  edm::Handle<std::vector<pat::Jet> > bottomCSVJets;
  evt.getByLabel("tightBottomCSVPFJets", bottomCSVJets);
  bottomCSVJetMultiplicity = bottomCSVJets.isValid() ? bottomCSVJets->size() : -1;

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
  
  edm::Handle< edm::View<reco::RecoCandidate> > leps;
  evt.getByLabel(leps_, leps);
  
  edm::Handle< edm::View<reco::RecoCandidate> > mets;
  evt.getByLabel(mets_, mets);
  
  edm::Handle<edm::View<PileupSummaryInfo> > PU_h;
  evt.getByLabel(PUSrc_, PU_h);
  
  if(PU_h.isValid()){
    for(edm::View<PileupSummaryInfo>::const_iterator iterPU = PU_h->begin(); iterPU != PU_h->end(); ++iterPU){  // vector size is 3
      int BX = iterPU->getBunchCrossing(); // -1: previous BX, 0: current BX,  1: next BX
      
      if      (BX == -1) nPU[0] = iterPU->getPU_NumInteractions();
      else if (BX ==  0) nPU[1] = iterPU->getPU_NumInteractions();
      else if (BX ==  1) nPU[2] = iterPU->getPU_NumInteractions();
    }
  }
  
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
  
  edm::Handle<double> mcWeightSrc_h;
  evt.getByLabel(mcWeightSrc_, mcWeightSrc_h);
  mcWeight = mcWeightSrc_h.isValid() ? *mcWeightSrc_h : -100.;
  
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
    for (int i = 0; i < 4; ++i) permutation[i] = jetLeptonCombinationCurrent[i];
    
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
      
      /* Permutation monitoring
      std::cout << "\nevent=" << evt.eventAuxiliary().event() << " hypoClassKey=" << hypoClassKey << " combi=" << h << " target=" << target << std::endl;
      std::cout << "jetLeptonCombinationGenMatch " << "jetLeptonCombinationCurrent[i] " << "jetLeptonCombinationCurrent2[i] " << "CSV " << std::endl;
      for (int i = 0; i < 4; ++i) {
        std::cout << jetLeptonCombinationGenMatch[i] << " " << jetLeptonCombinationCurrent[i] << " " << jetLeptonCombinationCurrent2[i] << " " << jets->at(i).bDiscriminator("combinedSecondaryVertexBJetTags") << std::endl;
      }
      //*/
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
    const reco::Candidate* hadW       = semiLepEvt->hadronicDecayW         (hypoClassKey, h);
    const reco::Candidate* hadB       = semiLepEvt->hadronicDecayB         (hypoClassKey, h);
    const reco::Candidate* hadQ       = semiLepEvt->hadronicDecayQuark     (hypoClassKey, h);
    const reco::Candidate* hadQBar    = semiLepEvt->hadronicDecayQuarkBar  (hypoClassKey, h);
    
    const reco::Candidate* lepTop      = semiLepEvt->leptonicDecayTop       (hypoClassKey, h);
    const reco::Candidate* lepW        = semiLepEvt->leptonicDecayW         (hypoClassKey, h);
    const reco::Candidate* lepB        = semiLepEvt->leptonicDecayB         (hypoClassKey, h);
    const reco::Candidate* lepton      = semiLepEvt->singleLepton           (hypoClassKey, h);
    const reco::Candidate* nu          = semiLepEvt->singleNeutrino         (hypoClassKey, h);
    
    hadQ_    = hadQ    ->p4();
    hadQBar_ = hadQBar ->p4();
    hadB_    = hadB    ->p4();
    lepB_    = lepB    ->p4();
    lepton_  = lepton  ->p4();
    nu_      = nu      ->p4();
    hadW_    = hadW    ->p4();
    lepW_    = lepW    ->p4();
    hadTop_  = hadTop  ->p4();
    lepTop_  = lepTop  ->p4();

    
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // get genParticles
    //////////////////////////////////////////////////////////////////////////////////////////////////

    const reco::Candidate* hadTopGen  = semiLepEvt->hadronicDecayTop();
    const reco::Candidate* hadBGen    = semiLepEvt->hadronicDecayB();
    const reco::Candidate* hadWGen    = semiLepEvt->hadronicDecayW();
    const reco::Candidate* hadQGen    = semiLepEvt->hadronicDecayQuark();
    const reco::Candidate* hadQBarGen = semiLepEvt->hadronicDecayQuarkBar();
    
    const reco::Candidate* lepTopGen  = semiLepEvt->leptonicDecayTop();
    const reco::Candidate* lepBGen    = semiLepEvt->leptonicDecayB();
    const reco::Candidate* lepWGen    = semiLepEvt->leptonicDecayW();
    const reco::Candidate* leptonGen  = semiLepEvt->singleLepton();
    const reco::Candidate* nuGen      = semiLepEvt->singleNeutrino();
    
    if (hadQGen && hadQBarGen && hadBGen && lepBGen && leptonGen && nuGen && hadWGen && lepWGen && hadTopGen && lepTopGen) {
      hadQGen_    = hadQGen    ->p4();
      hadQBarGen_ = hadQBarGen ->p4();
      hadBGen_    = hadBGen    ->p4();
      lepBGen_    = lepBGen    ->p4();
      leptonGen_  = leptonGen  ->p4();
      nuGen_      = nuGen      ->p4();
      hadWGen_    = hadWGen    ->p4();
      lepWGen_    = lepWGen    ->p4();
      hadTopGen_  = hadTopGen  ->p4();
      lepTopGen_  = lepTopGen  ->p4();
    }
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // get reconstructed particles (hypothesis might change particle momenta)
    //////////////////////////////////////////////////////////////////////////////////////////////////
    
    hadQRaw_    = jets->at(jetLeptonCombinationCurrent[0]).p4(); //correctedP4("L7Parton", "uds");
    hadQBarRaw_ = jets->at(jetLeptonCombinationCurrent[1]).p4(); //correctedP4("L7Parton", "uds");
    hadBRaw_    = jets->at(jetLeptonCombinationCurrent[2]).p4(); //correctedP4("L7Parton", "bottom");
    lepBRaw_    = jets->at(jetLeptonCombinationCurrent[3]).p4(); //correctedP4("L7Parton", "bottom");
    
    leptonRaw_  = leps->at(0).p4();
    nuRaw_      = mets->at(0).p4();
    
    hadWRaw_ = hadQRaw_ + hadQBarRaw_;
    lepWRaw_ = leptonRaw_ + nuRaw_;
    
    hadTopRaw_ = hadWRaw_ + hadBRaw_;
    lepTopRaw_ = lepWRaw_ + lepBRaw_;

    
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // fill eventTree with pt, eta and the masses of the reconstructed particles
    //////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Legacy
    
    hadQEta       = hadQ->eta();
    hadQRawPt     = hadQRaw_.Pt();
    hadQBarRawPt  = hadQBarRaw_.Pt();
    hadWRawMass   = hadWRaw_.M();
    hadBEta       = hadB->eta();
    hadBRawPt     = hadBRaw_.Pt();
    leptonPt      = lepton->pt();
    leptonEta     = lepton->eta();
    nuRawPt       = nuRaw_.Pt();
    lepWRawMass   = lepWRaw_.M();
    hadTopMass    = hadTop->mass();
    hadTopRawMass = hadTopRaw_.M();
    lepTopRawMass = lepTopRaw_.M();
    
    deltaRHadQHadQBar     = ROOT::Math::VectorUtil::DeltaR(hadQ->polarP4(), hadQBar->polarP4());
    deltaThetaHadQHadQBar = ROOT::Math::VectorUtil::Angle(hadQ->polarP4(), hadQBar->polarP4());
    deltaRHadWHadB        = ROOT::Math::VectorUtil::DeltaR(hadW->polarP4(), hadB->polarP4());
    deltaThetaHadWHadB    = ROOT::Math::VectorUtil::Angle(hadW->polarP4(), hadB->polarP4());
    if (hadWGen && hadBGen) genDeltaThetaHadWHadB = ROOT::Math::VectorUtil::Angle(hadWGen->polarP4(), hadBGen->polarP4());
    deltaRLepBLepton      = ROOT::Math::VectorUtil::DeltaR(lepton->polarP4(), lepB->polarP4());
    deltaThetaLepBLepton  = ROOT::Math::VectorUtil::Angle(lepton->polarP4(), lepB->polarP4());
    deltaRHadBLepB        = ROOT::Math::VectorUtil::DeltaR(hadB->polarP4(), lepB->polarP4());
    deltaThetaHadBLepB    = ROOT::Math::VectorUtil::Angle(hadB->polarP4(), lepB->polarP4());
    
    // Jet properties
    
    hadQBSSV   = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LightQ]).bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
    hadQBCSV   = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LightQ]).bDiscriminator("combinedSecondaryVertexBJetTags");
    hadQJC     = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LightQ]).jetCharge();
    hadQF      = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LightQ]).partonFlavour();
    
    hadQBarBSSV   = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LightQBar]).bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
    hadQBarBCSV   = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LightQBar]).bDiscriminator("combinedSecondaryVertexBJetTags");
    hadQBarJC     = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LightQBar]).jetCharge();
    hadQBarF      = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LightQBar]).partonFlavour();
    
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
  
    hadBBSSV   = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::HadB]).bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
    hadBBCSV   = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::HadB]).bDiscriminator("combinedSecondaryVertexBJetTags");
    hadBJC     = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::HadB]).jetCharge();
    hadBF      = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::HadB]).partonFlavour();
    
    leptonC    = leps->at(0).charge();
    leptonId   = (!leps_.label().compare("tightMuons"))?13:11;
    
    lepBBSSV   = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LepB]).bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
    lepBBCSV   = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LepB]).bDiscriminator("combinedSecondaryVertexBJetTags");
    lepBJC     = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LepB]).jetCharge();
    lepBF      = jets->at(jetLeptonCombinationCurrent[TtSemiLepEvtPartons::LepB]).partonFlavour();
    
    // Permutation properties
    
    mvaDisc    = semiLepEvt->mvaDisc(hMVA);
    fitChi2    = semiLepEvt->fitChi2(hKinFit);
    fitProb    = semiLepEvt->fitProb(hKinFit);
    hitFitChi2 = semiLepEvt->hitFitChi2(hHitFit);
    hitFitProb = semiLepEvt->hitFitProb(hHitFit);
    hitFitMT   = semiLepEvt->hitFitMT(hHitFit);
    hitFitSigMT= semiLepEvt->hitFitSigMT(hHitFit);
    
    const ROOT::Math::Boost CoMBoostHadTop(hadTop_.BoostToCM());
    const ROOT::Math::Boost CoMBoostLepTop(lepTop_.BoostToCM());
    
    amwtHadQ    = TMath::Max(0., 4*hadTop_.M()
                  *CoMBoostHadTop(hadQ_).E()*(hadTop_.M2()-hadB_.M2()
                    -2*hadTop_.M()*CoMBoostHadTop(hadQ_).E())
                  /(TMath::Power(hadTop_.M2()-hadB_.M2(),2)+hadW_.M2()*(hadTop_.M2()
                    -hadB_.M2())-TMath::Power(hadW_.M2(),2)));
    amwtHadQBar = TMath::Max(0., 4*hadTop_.M()
                  *CoMBoostHadTop(hadQBar_).E()*(hadTop_.M2()-hadB_.M2()
                    -2*hadTop_.M()*CoMBoostHadTop(hadQBar_).E())
                  /(TMath::Power(hadTop_.M2()-hadB_.M2(),2)+hadW_.M2()*(hadTop_.M2()
                    -hadB_.M2())-TMath::Power(hadW_.M2(),2)));
    amwtLepton  = TMath::Max(0., 4*lepTop_.M()
                  *CoMBoostLepTop(lepton_).E()*(lepTop_.M2()-lepB_.M2()
                    -2*lepTop_.M()*CoMBoostLepTop(lepton_).E())
                  /(TMath::Power(lepTop_.M2()-lepB_.M2(),2)+lepW_.M2()*(lepTop_.M2()
                    -lepB_.M2())-TMath::Power(lepW_.M2(),2)));
    
    // CR observables
    //const ROOT::Math::Boost CoMBoostHadW(hadW_.BoostToCM());
    
    inWEta = (hadQRaw_.Eta() + hadQBarRaw_.Eta())/2.;
    inWPhi = (hadQRaw_.Phi() + hadQBarRaw_.Phi())/2.;
    //std::cout << "abs(hadQRaw_.Phi() - hadQBarRaw_.Phi()): " << abs(hadQRaw_.Phi() - hadQBarRaw_.Phi()) << " TMath::Pi(): " << TMath::Pi() << std::endl;
    if (TMath::Abs(hadQRaw_.Phi() - hadQBarRaw_.Phi()) > TMath::Pi()) {
      if (inWPhi < 0.) {
        inWPhi = inWPhi + TMath::Pi();
      }
      else if (inWPhi > 0.) {
        inWPhi = inWPhi - TMath::Pi();
      }
    }
    
    hadWGeom_.SetPt(0.);
    hadWGeom_.SetEta(inWEta);
    hadWGeom_.SetPhi(inWPhi);
    hadWGeom_.SetM(0.);
    
    //std::cout << "hadQRaw_.Phi(): " << hadQRaw_.Phi() << " hadQBarRaw_.Phi(): " << hadQBarRaw_.Phi() << " inWPhi: " << inWPhi << std::endl;
    inWdR  = ROOT::Math::VectorUtil::DeltaR(hadQRaw_, hadQBarRaw_);
    
    //*
    inWJetMultiplicity = 0;
    inWContamination   = 0;
    if (   ROOT::Math::VectorUtil::DeltaR(hadWGeom_, hadBRaw_) < inWdR * 0.51
        || ROOT::Math::VectorUtil::DeltaR(hadWGeom_, lepBRaw_) < inWdR * 0.51) {
        inWContamination += 100;
    }
    for (int i = 0; i < allJetMultiplicity; ++i) {
      if (ROOT::Math::VectorUtil::DeltaR(hadWGeom_, allJets->at(i).p4()) < inWdR * 0.51) {
        //std::cout << "test" << std::endl;
        inWJetsPt[inWJetMultiplicity]   = allJets->at(i).pt();
        inWJetsEta[inWJetMultiplicity]  = allJets->at(i).eta();
        inWJetsPhi[inWJetMultiplicity]  = allJets->at(i).phi();
        inWJetsDRQ[inWJetMultiplicity]  = ROOT::Math::VectorUtil::DeltaR(hadQ_, allJets->at(i).p4());
        inWJetsDRQBar[inWJetMultiplicity]  = ROOT::Math::VectorUtil::DeltaR(hadQBar_, allJets->at(i).p4());
        inWJetsCHM[inWJetMultiplicity]  = allJets->at(i).chargedHadronMultiplicity();
        inWJetsNC[inWJetMultiplicity]   = allJets->at(i).nConstituents();
        if (allJets->at(i).pt() > 20) inWContamination += 10;
        if (allJets->at(i).pt() > 10) inWContamination +=  1;
        ++inWJetMultiplicity;
      }
    }
    //*/
    //std::cout << inWJetMultiplicity << std::endl;
    
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
  
  //*
  eventTree->Branch("hadTop.", "PxPyPzEVector", &hadTop_);
  eventTree->Branch("hadTopRaw.", "PxPyPzEVector", &hadTopRaw_);
  eventTree->Branch("hadTopGen.", "PxPyPzEVector", &hadTopGen_);
  
  eventTree->Branch("hadB.", "PxPyPzEVector", &hadB_);
  eventTree->Branch("hadBRaw.", "PxPyPzEVector", &hadBRaw_);
  eventTree->Branch("hadBGen.", "PxPyPzEVector", &hadBGen_);

  eventTree->Branch("hadW.", "PxPyPzEVector", &hadW_);
  eventTree->Branch("hadWRaw.", "PxPyPzEVector", &hadWRaw_);
  eventTree->Branch("hadWGen.", "PxPyPzEVector", &hadWGen_);
  
  eventTree->Branch("hadQ.", "PxPyPzEVector", &hadQ_);
  eventTree->Branch("hadQRaw.", "PxPyPzEVector", &hadQRaw_);
  eventTree->Branch("hadQGen.", "PxPyPzEVector", &hadQGen_);
  
  eventTree->Branch("hadQBar.", "PxPyPzEVector", &hadQBar_);
  eventTree->Branch("hadQBarRaw.", "PxPyPzEVector", &hadQBarRaw_);
  eventTree->Branch("hadQBarGen.", "PxPyPzEVector", &hadQBarGen_);
  
  eventTree->Branch("lepTop.", "PxPyPzEVector", &lepTop_);
  eventTree->Branch("lepTopRaw.", "PxPyPzEVector", &lepTopRaw_);
  eventTree->Branch("lepTopGen.", "PxPyPzEVector", &lepTopGen_);
  
  eventTree->Branch("lepB.", "PxPyPzEVector", &lepB_);
  eventTree->Branch("lepBRaw.", "PxPyPzEVector", &lepBRaw_);
  eventTree->Branch("lepBGen.", "PxPyPzEVector", &lepBGen_);

  eventTree->Branch("lepW.", "PxPyPzEVector", &lepW_);
  eventTree->Branch("lepWRaw.", "PxPyPzEVector", &lepWRaw_);
  eventTree->Branch("lepWGen.", "PxPyPzEVector", &lepWGen_);
  
  eventTree->Branch("lepton.", "PxPyPzEVector", &lepton_);
  eventTree->Branch("leptonRaw.", "PxPyPzEVector", &leptonRaw_);
  eventTree->Branch("leptonGen.", "PxPyPzEVector", &leptonGen_);
  
  eventTree->Branch("nu.", "PxPyPzEVector", &nu_);
  eventTree->Branch("nuRaw.", "PxPyPzEVector", &nuRaw_);
  eventTree->Branch("nuGen.", "PxPyPzEVector", &nuGen_);
  //*/
  
  eventTree->Branch("hadQEta", &hadQEta, "hadQEta/D");
  eventTree->Branch("hadQBSSV", &hadQBSSV, "hadQBSSV/D");
  eventTree->Branch("hadQBCSV", &hadQBCSV, "hadQBCSV/D");
  eventTree->Branch("hadQJC", &hadQJC, "hadQJC/D");
  eventTree->Branch("hadQF", &hadQF, "hadQF/D");
  eventTree->Branch("hadQRawPt", &hadQRawPt, "hadQRawPt/D");
  
  eventTree->Branch("hadQBarBSSV", &hadQBarBSSV, "hadQBarBSSV/D");
  eventTree->Branch("hadQBarBCSV", &hadQBarBCSV, "hadQBarBCSV/D");
  eventTree->Branch("hadQBarJC", &hadQBarJC, "hadQBarJC/D");
  eventTree->Branch("hadQBarF", &hadQBarF, "hadQBarF/D");
  eventTree->Branch("hadQBarRawPt", &hadQBarRawPt, "hadQBarRawPt/D");
  
  eventTree->Branch("hadWRawMass", &hadWRawMass, "hadWRawMass/D");
  
  eventTree->Branch("hadBEta", &hadBEta, "hadBEta/D");
  eventTree->Branch("hadBBSSV", &hadBBSSV, "hadBBSSV/D");
  eventTree->Branch("hadBBCSV", &hadBBCSV, "hadBBCSV/D");
  eventTree->Branch("hadBJC", &hadBJC, "hadBJC/D");
  eventTree->Branch("hadBF", &hadBF, "hadBF/D");
  
  eventTree->Branch("hadBRawPt", &hadBRawPt, "hadBRawPt/D");
  
  eventTree->Branch("leptonPt", &leptonPt, "leptonPt/D");
  eventTree->Branch("leptonEta", &leptonEta, "leptonEta/D");
  eventTree->Branch("leptonC", &leptonC, "leptonC/D");
  eventTree->Branch("leptonId", &leptonId, "leptonId/I");
  
  eventTree->Branch("nuRawPt", &nuRawPt, "nuRawPt/D");
  eventTree->Branch("lepWRawMass", &lepWRawMass, "lepWRawMass/D");
  
  eventTree->Branch("lepBBSSV", &lepBBSSV, "lepBBSSV/D");
  eventTree->Branch("lepBBCSV", &lepBBCSV, "lepBBCSV/D");
  eventTree->Branch("lepBJC", &lepBJC, "lepBJC/D");
  eventTree->Branch("lepBF", &lepBF, "lepBF/D");
  
  eventTree->Branch("hadTopMass", &hadTopMass, "hadTopMass/D");
  eventTree->Branch("hadTopRawMass", &hadTopRawMass, "hadTopRawMass/D");  
  eventTree->Branch("lepTopRawMass", &lepTopRawMass, "lepTopRawMass/D");
  
  eventTree->Branch("deltaRHadQHadQBar", &deltaRHadQHadQBar, "deltaRHadQHadQBar/D");
  eventTree->Branch("deltaThetaHadQHadQBar", &deltaThetaHadQHadQBar, "deltaThetaHadQHadQBar/D");
  eventTree->Branch("deltaRHadWHadB", &deltaRHadWHadB, "deltaRHadWHadB/D");
  eventTree->Branch("deltaThetaHadWHadB", &deltaThetaHadWHadB, "deltaThetaHadWHadB/D");
  eventTree->Branch("genDeltaThetaHadWHadB", &genDeltaThetaHadWHadB, "genDeltaThetaHadWHadB/D");
  eventTree->Branch("deltaRLepBLepton", &deltaRLepBLepton, "deltaRLepBLepton/D");
  eventTree->Branch("deltaThetaLepBLepton", &deltaThetaLepBLepton, "deltaThetaLepBLepton/D");
  eventTree->Branch("deltaRHadBLepB", &deltaRHadBLepB, "deltaRHadBLepB/D");
  eventTree->Branch("deltaThetaHadBLepB", &deltaThetaHadBLepB, "deltaThetaHadBLepB/D");
  
  eventTree->Branch("jetMultiplicity", &jetMultiplicity, "jetMultiplicity/I");
  eventTree->Branch("jetsPt", &jetsPt, "jetsPt[6]/D");
  eventTree->Branch("jetsEta", &jetsEta, "jetsEta[6]/D");
  eventTree->Branch("jetsPhi", &jetsPhi, "jetsPhi[6]/D");
  eventTree->Branch("jetsBCSV", &jetsBCSV, "jetsBCSV[6]/D");
  eventTree->Branch("permutation", &permutation, "permutation[4]/I");
  eventTree->Branch("noPtEtaJetMultiplicity", &noPtEtaJetMultiplicity, "noPtEtaJetMultiplicity/I");
  eventTree->Branch("bottomSSVJetMultiplicity", &bottomSSVJetMultiplicity, "bottomSSVJetMultiplicity/I");
  eventTree->Branch("bottomCSVJetMultiplicity", &bottomCSVJetMultiplicity, "bottomCSVJetMultiplicity/I");
  
  eventTree->Branch("allJetMultiplicity", &allJetMultiplicity, "allJetMultiplicity/I");
  eventTree->Branch("inWJetMultiplicity", &inWJetMultiplicity, "inWJetMultiplicity/I");
  eventTree->Branch("inWContamination", &inWContamination, "inWContamination/I");
  eventTree->Branch("inWJetsPt", &inWJetsPt, "inWJetsPt[inWJetMultiplicity]/D");
  eventTree->Branch("inWJetsEta", &inWJetsEta, "inWJetsEta[inWJetMultiplicity]/D");
  eventTree->Branch("inWJetsPhi", &inWJetsPhi, "inWJetsPhi[inWJetMultiplicity]/D");
  eventTree->Branch("inWJetsDRQ", &inWJetsDRQ, "inWJetsDRQ[inWJetMultiplicity]/D");
  eventTree->Branch("inWJetsDRQBar", &inWJetsDRQBar, "inWJetsDRQBar[inWJetMultiplicity]/D");
  eventTree->Branch("inWJetsCHM", &inWJetsCHM, "inWJetsCHM[inWJetMultiplicity]/D");
  eventTree->Branch("inWJetsNC", &inWJetsNC, "inWJetsNC[inWJetMultiplicity]/D");
  
  eventTree->Branch("genMatchDr", &genMatchDr, "genMatchDr/D");
  eventTree->Branch("mvaDisc", &mvaDisc, "mvaDisc/D");
  eventTree->Branch("fitChi2", &fitChi2, "fitChi2/D");
  eventTree->Branch("fitProb", &fitProb, "fitProb/D");
  eventTree->Branch("hitFitChi2", &hitFitChi2, "hitFitChi2/D");
  eventTree->Branch("hitFitProb", &hitFitProb, "hitFitProb/D");
  eventTree->Branch("hitFitMT", &hitFitMT, "hitFitMT/D");
  eventTree->Branch("hitFitSigMT", &hitFitSigMT, "hitFitSigMT/D");
  
  eventTree->Branch("amwtHadQ", &amwtHadQ, "amwtHadQ/D");
  eventTree->Branch("amwtHadQBar", &amwtHadQBar, "amwtHadQBar/D");
  eventTree->Branch("amwtLepton", &amwtLepton, "amwtLepton/D");
  
  eventTree->Branch("inWEta", &inWEta, "inWEta/D");
  eventTree->Branch("inWPhi", &inWPhi, "inWPhi/D");
  eventTree->Branch("inWdR", &inWdR, "inWdR/D");
  
  eventTree->Branch("nPU", &nPU, "nPU[3]/I");
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
  eventTree->Branch("mcWeight", &mcWeight, "mcWeight/D");
  eventTree->Branch("MCWeight", &MCWeight, "MCWeight/D");
  
  if(savePDFWeights_) {
    LHAPDF::initPDFSet(1, "cteq66.LHgrid");
    eventTree->Branch("pdfWeights", &pdfWeights, "pdfWeights[44]/D");
  }
  
  eventTree->Branch("target", &target, "target/I");
}

void
EventHypothesisAnalyzer::endJob() 
{
}

EventHypothesisAnalyzer::~EventHypothesisAnalyzer() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(EventHypothesisAnalyzer);
