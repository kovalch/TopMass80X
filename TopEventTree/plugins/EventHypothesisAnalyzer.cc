#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


#include "AnalysisDataFormats/TopObjects/interface/TtSemiLepEvtPartons.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "LHAPDF/LHAPDF.h"

#include "TMath.h"

#include "TopMass/TopEventTree/plugins/EventHypothesisAnalyzer.h"

EventHypothesisAnalyzer::EventHypothesisAnalyzer(const edm::ParameterSet& cfg):
  eventTree(0),
  ttEvent_  (cfg.getParameter<edm::InputTag>("ttEvent")),
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
  savePDFWeights_(cfg.getParameter<bool>("savePDFWeights" )),
  
  kJetMAX(cfg.getParameter<int>("maxNJets" )),

  treeToAppend_(cfg.getParameter<std::string>("treeToAppend"))
{
}

void
EventHypothesisAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& setup)
{
  //////////////////////////////////////////////////////////////////////////////
  // get a handle for the TtSemiLeptonicEvent and a key to the hypothesis
  ////////////////////////////////////////////////////////////////////////////
  
  edm::Handle<TtSemiLeptonicEvent> ttEvent;
  evt.getByLabel(ttEvent_, ttEvent);
  
  edm::Handle<int> hypoClassKeyHandle;
  evt.getByLabel(hypoClassKey_, hypoClassKeyHandle);
  TtSemiLeptonicEvent::HypoClassKey& hypoClassKey = (TtSemiLeptonicEvent::HypoClassKey&) *hypoClassKeyHandle;
  TtSemiLeptonicEvent::HypoClassKey hypoClassKeyGenMatch = (TtSemiLeptonicEvent::HypoClassKey) 3;
  TtSemiLeptonicEvent::HypoClassKey hypoClassKeyMVA      = (TtSemiLeptonicEvent::HypoClassKey) 4;
  TtSemiLeptonicEvent::HypoClassKey hypoClassKeyKinFit   = (TtSemiLeptonicEvent::HypoClassKey) 5;
  TtSemiLeptonicEvent::HypoClassKey hypoClassKeyHitFit   = (TtSemiLeptonicEvent::HypoClassKey) 8;
  
  //////////////////////////////////////////////////////////////////////////////
  // INIT TopEvent
  ////////////////////////////////////////////////////////////////////////////
  
  top->init();
  
  top->run       = evt.eventAuxiliary().run();
  top->lumiBlock = evt.eventAuxiliary().luminosityBlock();
  top->event     = evt.eventAuxiliary().event();
  
  //////////////////////////////////////////////////////////////////////////////
  // JETS
  ////////////////////////////////////////////////////////////////////////////
  
  edm::Handle<std::vector<pat::Jet> > jets;
  evt.getByLabel(jets_, jets);
  
  for(std::vector< pat::Jet >::const_iterator ijet = jets->begin(); ijet != jets->end(); ++ijet) {
    
    top->jet.push_back(TLorentzVector(ijet->px(), ijet->py(), ijet->pz(), ijet->energy()));
    
    top->jetCharge.push_back(ijet->jetCharge());
    top->jetFlavour.push_back(ijet->partonFlavour());
    top->jetCSV.push_back(ijet->bDiscriminator("combinedSecondaryVertexBJetTags"));
  }
  

  
  //////////////////////////////////////////////////////////////////////////////
  // GENPARTON
  // TTBar, HadTop, HadW, LepTop, LepW,
  // HadB, LightQ, LightQBar, LepB, Lepton, Neutrino
  ////////////////////////////////////////////////////////////////////////////
  
  if (ttEvent->hadronicDecayTop()) {
    top->genparton.push_back(TLorentzVector(
      ttEvent->hadronicDecayTop()->px()     + ttEvent->leptonicDecayTop()->px(),
      ttEvent->hadronicDecayTop()->py()     + ttEvent->leptonicDecayTop()->py(),
      ttEvent->hadronicDecayTop()->pz()     + ttEvent->leptonicDecayTop()->pz(),
      ttEvent->hadronicDecayTop()->energy() + ttEvent->leptonicDecayTop()->energy()));
    
    top->genparton.push_back(TLorentzVector(
      ttEvent->hadronicDecayTop()->px(),      ttEvent->hadronicDecayTop()->py(),
      ttEvent->hadronicDecayTop()->pz(),      ttEvent->hadronicDecayTop()->energy()));
    top->genparton.push_back(TLorentzVector(
      ttEvent->leptonicDecayTop()->px(),      ttEvent->leptonicDecayTop()->py(),
      ttEvent->leptonicDecayTop()->pz(),      ttEvent->leptonicDecayTop()->energy()));
    
    top->genparton.push_back(TLorentzVector(
      ttEvent->hadronicDecayW()->px(),        ttEvent->hadronicDecayW()->py(),
      ttEvent->hadronicDecayW()->pz(),        ttEvent->hadronicDecayW()->energy()));
    top->genparton.push_back(TLorentzVector(
      ttEvent->leptonicDecayW()->px(),        ttEvent->leptonicDecayW()->py(),
      ttEvent->leptonicDecayW()->pz(),        ttEvent->leptonicDecayW()->energy()));
    
    top->genparton.push_back(TLorentzVector(
      ttEvent->hadronicDecayB()->px(),        ttEvent->hadronicDecayB()->py(),
      ttEvent->hadronicDecayB()->pz(),        ttEvent->hadronicDecayB()->energy()));
    top->genparton.push_back(TLorentzVector(
      ttEvent->hadronicDecayQuark()->px(),    ttEvent->hadronicDecayQuark()->py(),
      ttEvent->hadronicDecayQuark()->pz(),    ttEvent->hadronicDecayQuark()->energy()));
    top->genparton.push_back(TLorentzVector(
      ttEvent->hadronicDecayQuarkBar()->px(), ttEvent->hadronicDecayQuarkBar()->py(),
      ttEvent->hadronicDecayQuarkBar()->pz(), ttEvent->hadronicDecayQuarkBar()->energy()));
    top->genparton.push_back(TLorentzVector(
      ttEvent->leptonicDecayB()->px(),        ttEvent->leptonicDecayB()->py(),
      ttEvent->leptonicDecayB()->pz(),        ttEvent->leptonicDecayB()->energy()));
    top->genparton.push_back(TLorentzVector(
      ttEvent->singleLepton()->px(),          ttEvent->singleLepton()->py(),
      ttEvent->singleLepton()->pz(),          ttEvent->singleLepton()->energy()));
    top->genparton.push_back(TLorentzVector(
      ttEvent->singleNeutrino()->px(),        ttEvent->singleNeutrino()->py(),
      ttEvent->singleNeutrino()->pz(),        ttEvent->singleNeutrino()->energy()));
  }
  
  if (ttEvent->isHypoValid(hypoClassKeyGenMatch)) {
    for (int p = 0; p < 11; ++p) {
      switch(p) {
        case HadB:
          top->genpartonJetIdx.push_back(ttEvent->jetLeptonCombination(hypoClassKeyGenMatch)[TtSemiLepEvtPartons::HadB]);
          break;
        case LightQ:
          top->genpartonJetIdx.push_back(ttEvent->jetLeptonCombination(hypoClassKeyGenMatch)[TtSemiLepEvtPartons::LightQ]);
          break;
        case LightQBar:
          top->genpartonJetIdx.push_back(ttEvent->jetLeptonCombination(hypoClassKeyGenMatch)[TtSemiLepEvtPartons::LightQBar]);
          break;
        case LepB:
          top->genpartonJetIdx.push_back(ttEvent->jetLeptonCombination(hypoClassKeyGenMatch)[TtSemiLepEvtPartons::LepB]);
          break;
        default: top->genpartonJetIdx.push_back(-1);
      }
    }
  }
  
  //////////////////////////////////////////////////////////////////////////////
  // FIT AND RECO
  ////////////////////////////////////////////////////////////////////////////
  
  for (unsigned int h = 0; h < ttEvent->numberOfAvailableHypos(hypoClassKey); ++h) {
    if (!ttEvent->isHypoValid(hypoClassKey, h)) break;
    std::vector<TLorentzVector> hfit = getPartons(ttEvent, hypoClassKey, h);
    top->fit.push_back(hfit);
    
    // TODO switch kinfit/hitfit
    top->fitProb.push_back(ttEvent->hitFitProb(h));
    top->fitChi2.push_back(ttEvent->hitFitChi2(h));
    
    std::vector<TLorentzVector> hreco = getPartons(ttEvent, hypoClassKeyMVA, 
      ttEvent->correspondingHypo(hypoClassKey, h, hypoClassKeyMVA));
    top->reco.push_back(hreco);
    
    std::vector<int> hrecoJetIdx;
    for (int p = 0; p < 11; ++p) {
      switch(p) {
        case HadB:
          hrecoJetIdx.push_back(ttEvent->jetLeptonCombination(hypoClassKey, h)[TtSemiLepEvtPartons::HadB]);
          break;
        case LightQ:
          hrecoJetIdx.push_back(ttEvent->jetLeptonCombination(hypoClassKey, h)[TtSemiLepEvtPartons::LightQ]);
          break;
        case LightQBar:
          hrecoJetIdx.push_back(ttEvent->jetLeptonCombination(hypoClassKey, h)[TtSemiLepEvtPartons::LightQBar]);
          break;
        case LepB:
          hrecoJetIdx.push_back(ttEvent->jetLeptonCombination(hypoClassKey, h)[TtSemiLepEvtPartons::LepB]);
          break;
        default: hrecoJetIdx.push_back(-1);
      }
    }
    top->recoJetIdx.push_back(hrecoJetIdx);
  }
  
  // Fetch reconstructed permutations without fit solution
  
  for (unsigned int h = 0; h < ttEvent->numberOfAvailableHypos(hypoClassKeyMVA); ++h) {
    if (!ttEvent->isHypoValid(hypoClassKeyMVA, h)) break;
    if ( ttEvent->correspondingHypo(hypoClassKeyMVA, h, hypoClassKey) > -1) continue;
    
    std::vector<TLorentzVector> hreco = getPartons(ttEvent, hypoClassKeyMVA, h);
    top->reco.push_back(hreco);
    
    std::vector<int> hrecoJetIdx;
    for (int p = 0; p < 11; ++p) {
      switch(p) {
        case HadB:
          hrecoJetIdx.push_back(ttEvent->jetLeptonCombination(hypoClassKeyMVA, h)[TtSemiLepEvtPartons::HadB]);
          break;
        case LightQ:
          hrecoJetIdx.push_back(ttEvent->jetLeptonCombination(hypoClassKeyMVA, h)[TtSemiLepEvtPartons::LightQ]);
          break;
        case LightQBar:
          hrecoJetIdx.push_back(ttEvent->jetLeptonCombination(hypoClassKeyMVA, h)[TtSemiLepEvtPartons::LightQBar]);
          break;
        case LepB:
          hrecoJetIdx.push_back(ttEvent->jetLeptonCombination(hypoClassKeyMVA, h)[TtSemiLepEvtPartons::LepB]);
          break;
        default: hrecoJetIdx.push_back(-1);
      }
    }
    top->recoJetIdx.push_back(hrecoJetIdx);
  }
  
  eventTree->Fill();
}

void 
EventHypothesisAnalyzer::beginJob()
{
  edm::Service<TFileService> fs;
  if( !fs ) throw edm::Exception( edm::errors::Configuration, "TFile Service is not registered in cfg file" );
  
  top = new TopEvent();
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // get tree
  //////////////////////////////////////////////////////////////////////////////////////////////////
  if(treeToAppend_.size()) eventTree = (TTree*)fs->file().Get(treeToAppend_.c_str());

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // create tree in case no tree was found
  //////////////////////////////////////////////////////////////////////////////////////////////////
  if(!eventTree) eventTree = fs->make<TTree>("eventTree", "Tree for UHH top-quark analysis\nParticles are in order {TTBar, HadTop, LepTop, HadW, LepW, HadB, LightQ, LightQBar, LepB, Lepton, Neutrino}");

  eventTree->Branch("top", top); //, 32000, -1);
}

std::vector<TLorentzVector> EventHypothesisAnalyzer::getPartons(edm::Handle<TtSemiLeptonicEvent> ttEvent, TtSemiLeptonicEvent::HypoClassKey hypoClassKey, unsigned int h) {
  // TODO validate hypo, catch h=-1
  std::vector<TLorentzVector> parton;

  parton.push_back(TLorentzVector(
    ttEvent->hadronicDecayTop(hypoClassKey, h)->px()     + ttEvent->leptonicDecayTop(hypoClassKey, h)->px(),
    ttEvent->hadronicDecayTop(hypoClassKey, h)->py()     + ttEvent->leptonicDecayTop(hypoClassKey, h)->py(),
    ttEvent->hadronicDecayTop(hypoClassKey, h)->pz()     + ttEvent->leptonicDecayTop(hypoClassKey, h)->pz(),
    ttEvent->hadronicDecayTop(hypoClassKey, h)->energy() + ttEvent->leptonicDecayTop(hypoClassKey, h)->energy()));
  
  parton.push_back(TLorentzVector(
    ttEvent->hadronicDecayTop(hypoClassKey, h)->px(),      ttEvent->hadronicDecayTop(hypoClassKey, h)->py(),
    ttEvent->hadronicDecayTop(hypoClassKey, h)->pz(),      ttEvent->hadronicDecayTop(hypoClassKey, h)->energy()));
  parton.push_back(TLorentzVector(
    ttEvent->leptonicDecayTop(hypoClassKey, h)->px(),      ttEvent->leptonicDecayTop(hypoClassKey, h)->py(),
    ttEvent->leptonicDecayTop(hypoClassKey, h)->pz(),      ttEvent->leptonicDecayTop(hypoClassKey, h)->energy()));
  
  parton.push_back(TLorentzVector(
    ttEvent->hadronicDecayW(hypoClassKey, h)->px(),        ttEvent->hadronicDecayW(hypoClassKey, h)->py(),
    ttEvent->hadronicDecayW(hypoClassKey, h)->pz(),        ttEvent->hadronicDecayW(hypoClassKey, h)->energy()));
  parton.push_back(TLorentzVector(
    ttEvent->leptonicDecayW(hypoClassKey, h)->px(),        ttEvent->leptonicDecayW(hypoClassKey, h)->py(),
    ttEvent->leptonicDecayW(hypoClassKey, h)->pz(),        ttEvent->leptonicDecayW(hypoClassKey, h)->energy()));
  
  parton.push_back(TLorentzVector(
    ttEvent->hadronicDecayB(hypoClassKey, h)->px(),        ttEvent->hadronicDecayB(hypoClassKey, h)->py(),
    ttEvent->hadronicDecayB(hypoClassKey, h)->pz(),        ttEvent->hadronicDecayB(hypoClassKey, h)->energy()));
  parton.push_back(TLorentzVector(
    ttEvent->hadronicDecayQuark(hypoClassKey, h)->px(),    ttEvent->hadronicDecayQuark(hypoClassKey, h)->py(),
    ttEvent->hadronicDecayQuark(hypoClassKey, h)->pz(),    ttEvent->hadronicDecayQuark(hypoClassKey, h)->energy()));
  parton.push_back(TLorentzVector(
    ttEvent->hadronicDecayQuarkBar(hypoClassKey, h)->px(), ttEvent->hadronicDecayQuarkBar(hypoClassKey, h)->py(),
    ttEvent->hadronicDecayQuarkBar(hypoClassKey, h)->pz(), ttEvent->hadronicDecayQuarkBar(hypoClassKey, h)->energy()));
  parton.push_back(TLorentzVector(
    ttEvent->leptonicDecayB(hypoClassKey, h)->px(),        ttEvent->leptonicDecayB(hypoClassKey, h)->py(),
    ttEvent->leptonicDecayB(hypoClassKey, h)->pz(),        ttEvent->leptonicDecayB(hypoClassKey, h)->energy()));
  parton.push_back(TLorentzVector(
    ttEvent->singleLepton(hypoClassKey, h)->px(),          ttEvent->singleLepton(hypoClassKey, h)->py(),
    ttEvent->singleLepton(hypoClassKey, h)->pz(),          ttEvent->singleLepton(hypoClassKey, h)->energy()));
  parton.push_back(TLorentzVector(
    ttEvent->singleNeutrino(hypoClassKey, h)->px(),        ttEvent->singleNeutrino(hypoClassKey, h)->py(),
    ttEvent->singleNeutrino(hypoClassKey, h)->pz(),        ttEvent->singleNeutrino(hypoClassKey, h)->energy()));
  
  return parton;
}

void
EventHypothesisAnalyzer::endJob() 
{
}

EventHypothesisAnalyzer::~EventHypothesisAnalyzer() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(EventHypothesisAnalyzer);
