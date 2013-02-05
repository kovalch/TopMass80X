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

#include "TopMass/TopEventTree/interface/TreeRegistryService.h"

#include "AnalysisDataFormats/TopObjects/interface/TtSemiLepEvtPartons.h"
#include "AnalysisDataFormats/TopObjects/interface/TtFullHadEvtPartons.h"
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
savePDFWeights_(cfg.getParameter<bool>("savePDFWeights")),

kJetMAX_(cfg.getParameter<int>("maxNJets")),
kMAXCombo_(cfg.getParameter<int>("maxCombo"))
{
}

void
EventHypothesisAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& setup)
{
  //////////////////////////////////////////////////////////////////////////////
  // INIT TopEvent
  ////////////////////////////////////////////////////////////////////////////

  top->init();

  top->run       = evt.eventAuxiliary().run();
  top->lumiBlock = evt.eventAuxiliary().luminosityBlock();
  top->event     = evt.eventAuxiliary().event();

  //////////////////////////////////////////////////////////////////////////////
  // get a handle for the TtSemiLeptonicEvent and a key to the hypothesis
  ////////////////////////////////////////////////////////////////////////////

  edm::Handle<TtSemiLeptonicEvent> hSemiLepTtEvent;
  edm::Handle<TtFullHadronicEvent> hFullHadTtEvent;
  const TtSemiLeptonicEvent *semiLepTtEvent = 0;
  const TtFullHadronicEvent *fullHadTtEvent = 0;
  const TtEvent *ttEvent = 0;

  if(ttEvent_.label().find("SemiLep")!=std::string::npos){
    evt.getByLabel(ttEvent_, hSemiLepTtEvent);
    semiLepTtEvent = hSemiLepTtEvent.product();
    ttEvent = semiLepTtEvent;
  }
  else if(ttEvent_.label().find("FullHad")!=std::string::npos){
    evt.getByLabel(ttEvent_, hFullHadTtEvent);
    fullHadTtEvent = hFullHadTtEvent.product();
    ttEvent = fullHadTtEvent;
  }
  else{
    std::cout << "The given 'ttEvent' label (" << ttEvent_.label() << ") is not allowed.\n"
        << "It has to contain either 'SemiLep' or 'FullHad'!" << std::endl;
  }

  edm::Handle<int> hypoClassKeyHandle;
  evt.getByLabel(hypoClassKey_, hypoClassKeyHandle);

  TtEvent::HypoClassKey& hypoClassKey = (TtEvent::HypoClassKey&) *hypoClassKeyHandle;

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


  ////////////////////////////////////////////////////////////////////////////////
  // GENPARTON SEMI-LEP
  // TTBar, HadTop, LepTop, HadW, LepW,
  // HadB, LightQ, LightQBar, LepB, Lepton, Neutrino
  ////////////////////////////////////////////////////////////////////////////

  if(semiLepTtEvent){
    if (semiLepTtEvent->hadronicDecayTop()) {
      top->genparton.push_back(TLorentzVector(
          semiLepTtEvent->hadronicDecayTop()->px()     + semiLepTtEvent->leptonicDecayTop()->px(),
          semiLepTtEvent->hadronicDecayTop()->py()     + semiLepTtEvent->leptonicDecayTop()->py(),
          semiLepTtEvent->hadronicDecayTop()->pz()     + semiLepTtEvent->leptonicDecayTop()->pz(),
          semiLepTtEvent->hadronicDecayTop()->energy() + semiLepTtEvent->leptonicDecayTop()->energy()));

      top->genparton.push_back(TLorentzVector(
          semiLepTtEvent->hadronicDecayTop()->px(), semiLepTtEvent->hadronicDecayTop()->py(),
          semiLepTtEvent->hadronicDecayTop()->pz(), semiLepTtEvent->hadronicDecayTop()->energy()));
      top->genparton.push_back(TLorentzVector(
          semiLepTtEvent->leptonicDecayTop()->px(), semiLepTtEvent->leptonicDecayTop()->py(),
          semiLepTtEvent->leptonicDecayTop()->pz(), semiLepTtEvent->leptonicDecayTop()->energy()));

      top->genparton.push_back(TLorentzVector(
          semiLepTtEvent->hadronicDecayW()->px(), semiLepTtEvent->hadronicDecayW()->py(),
          semiLepTtEvent->hadronicDecayW()->pz(), semiLepTtEvent->hadronicDecayW()->energy()));
      top->genparton.push_back(TLorentzVector(
          semiLepTtEvent->leptonicDecayW()->px(), semiLepTtEvent->leptonicDecayW()->py(),
          semiLepTtEvent->leptonicDecayW()->pz(), semiLepTtEvent->leptonicDecayW()->energy()));

      top->genparton.push_back(TLorentzVector(
          semiLepTtEvent->hadronicDecayB()->px(), semiLepTtEvent->hadronicDecayB()->py(),
          semiLepTtEvent->hadronicDecayB()->pz(), semiLepTtEvent->hadronicDecayB()->energy()));
      top->genparton.push_back(TLorentzVector(
          semiLepTtEvent->hadronicDecayQuark()->px(), semiLepTtEvent->hadronicDecayQuark()->py(),
          semiLepTtEvent->hadronicDecayQuark()->pz(), semiLepTtEvent->hadronicDecayQuark()->energy()));
      top->genparton.push_back(TLorentzVector(
          semiLepTtEvent->hadronicDecayQuarkBar()->px(), semiLepTtEvent->hadronicDecayQuarkBar()->py(),
          semiLepTtEvent->hadronicDecayQuarkBar()->pz(), semiLepTtEvent->hadronicDecayQuarkBar()->energy()));
      top->genparton.push_back(TLorentzVector(
          semiLepTtEvent->leptonicDecayB()->px(), semiLepTtEvent->leptonicDecayB()->py(),
          semiLepTtEvent->leptonicDecayB()->pz(), semiLepTtEvent->leptonicDecayB()->energy()));
      top->genparton.push_back(TLorentzVector(
          semiLepTtEvent->singleLepton()->px(), semiLepTtEvent->singleLepton()->py(),
          semiLepTtEvent->singleLepton()->pz(), semiLepTtEvent->singleLepton()->energy()));
      top->genparton.push_back(TLorentzVector(
          semiLepTtEvent->singleNeutrino()->px(), semiLepTtEvent->singleNeutrino()->py(),
          semiLepTtEvent->singleNeutrino()->pz(), semiLepTtEvent->singleNeutrino()->energy()));
    }

    if (semiLepTtEvent->isHypoValid(TtEvent::kGenMatch)) {
      for (int p = 0; p < 11; ++p) {
        switch(p) {
        case HadB:
          top->genpartonJetIdx.push_back(semiLepTtEvent->jetLeptonCombination(TtEvent::kGenMatch)[TtSemiLepEvtPartons::HadB]);
          break;
        case LightQ:
          top->genpartonJetIdx.push_back(semiLepTtEvent->jetLeptonCombination(TtEvent::kGenMatch)[TtSemiLepEvtPartons::LightQ]);
          break;
        case LightQBar:
          top->genpartonJetIdx.push_back(semiLepTtEvent->jetLeptonCombination(TtEvent::kGenMatch)[TtSemiLepEvtPartons::LightQBar]);
          break;
        case LepB:
          top->genpartonJetIdx.push_back(semiLepTtEvent->jetLeptonCombination(TtEvent::kGenMatch)[TtSemiLepEvtPartons::LepB]);
          break;
        default: top->genpartonJetIdx.push_back(-1);
        }
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////////
  // GENPARTON FULL-HAD
  // TTBar, Top1, Top2, W1, W2,
  // B1, LightQ1, LightQBar1, B2, LightQ2, LightQBar2
  ////////////////////////////////////////////////////////////////////////////

  else if(fullHadTtEvent){
    if (fullHadTtEvent->top()) {
      top->genparton.push_back(TLorentzVector(
          fullHadTtEvent->top()->px()     + fullHadTtEvent->topBar()->px(),
          fullHadTtEvent->top()->py()     + fullHadTtEvent->topBar()->py(),
          fullHadTtEvent->top()->pz()     + fullHadTtEvent->topBar()->pz(),
          fullHadTtEvent->top()->energy() + fullHadTtEvent->topBar()->energy()));

      top->genparton.push_back(TLorentzVector(
          fullHadTtEvent->top()->px(), fullHadTtEvent->top()->py(),
          fullHadTtEvent->top()->pz(), fullHadTtEvent->top()->energy()));
      top->genparton.push_back(TLorentzVector(
          fullHadTtEvent->topBar()->px(), fullHadTtEvent->topBar()->py(),
          fullHadTtEvent->topBar()->pz(), fullHadTtEvent->topBar()->energy()));

      top->genparton.push_back(TLorentzVector(
          fullHadTtEvent->wPlus()->px(), fullHadTtEvent->wPlus()->py(),
          fullHadTtEvent->wPlus()->pz(), fullHadTtEvent->wPlus()->energy()));
      top->genparton.push_back(TLorentzVector(
          fullHadTtEvent->wMinus()->px(), fullHadTtEvent->wMinus()->py(),
          fullHadTtEvent->wMinus()->pz(), fullHadTtEvent->wMinus()->energy()));

      top->genparton.push_back(TLorentzVector(
          fullHadTtEvent->b()->px(), fullHadTtEvent->b()->py(),
          fullHadTtEvent->b()->pz(), fullHadTtEvent->b()->energy()));
      top->genparton.push_back(TLorentzVector(
          fullHadTtEvent->lightQ()->px(), fullHadTtEvent->lightQ()->py(),
          fullHadTtEvent->lightQ()->pz(), fullHadTtEvent->lightQ()->energy()));
      top->genparton.push_back(TLorentzVector(
          fullHadTtEvent->lightQBar()->px(), fullHadTtEvent->lightQBar()->py(),
          fullHadTtEvent->lightQBar()->pz(), fullHadTtEvent->lightQBar()->energy()));
      top->genparton.push_back(TLorentzVector(
          fullHadTtEvent->bBar()->px(), fullHadTtEvent->bBar()->py(),
          fullHadTtEvent->bBar()->pz(), fullHadTtEvent->bBar()->energy()));
      top->genparton.push_back(TLorentzVector(
          fullHadTtEvent->lightP()->px(), fullHadTtEvent->lightP()->py(),
          fullHadTtEvent->lightP()->pz(), fullHadTtEvent->lightP()->energy()));
      top->genparton.push_back(TLorentzVector(
          fullHadTtEvent->lightPBar()->px(), fullHadTtEvent->lightPBar()->py(),
          fullHadTtEvent->lightPBar()->pz(), fullHadTtEvent->lightPBar()->energy()));
    }

    if (fullHadTtEvent->isHypoValid(TtEvent::kGenMatch)) {
      for (int p = 0; p < 11; ++p) {
        switch(p) {
        case B1:
          top->genpartonJetIdx.push_back(fullHadTtEvent->jetLeptonCombination(TtEvent::kGenMatch)[TtFullHadEvtPartons::B]);
          break;
        case LightQ1:
          top->genpartonJetIdx.push_back(fullHadTtEvent->jetLeptonCombination(TtEvent::kGenMatch)[TtFullHadEvtPartons::LightQ]);
          break;
        case LightQBar1:
          top->genpartonJetIdx.push_back(fullHadTtEvent->jetLeptonCombination(TtEvent::kGenMatch)[TtFullHadEvtPartons::LightQBar]);
          break;
        case B2:
          top->genpartonJetIdx.push_back(fullHadTtEvent->jetLeptonCombination(TtEvent::kGenMatch)[TtFullHadEvtPartons::BBar]);
          break;
        case LightQ2:
          top->genpartonJetIdx.push_back(fullHadTtEvent->jetLeptonCombination(TtEvent::kGenMatch)[TtFullHadEvtPartons::LightP]);
          break;
        case LightQBar2:
          top->genpartonJetIdx.push_back(fullHadTtEvent->jetLeptonCombination(TtEvent::kGenMatch)[TtFullHadEvtPartons::LightPBar]);
          break;
        default: top->genpartonJetIdx.push_back(-1);
        }
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // FIT AND RECO
  ////////////////////////////////////////////////////////////////////////////

  for (unsigned int h = 0, maxHypos = std::min(ttEvent->numberOfAvailableHypos(hypoClassKey), kMAXCombo_); h < maxHypos; ++h) {
    if (!ttEvent->isHypoValid(hypoClassKey, h)) break;
    std::vector<TLorentzVector> hfit;
    if     (semiLepTtEvent) hfit = getPartons(semiLepTtEvent, hypoClassKey, h);
    else if(fullHadTtEvent) hfit = getPartons(fullHadTtEvent, hypoClassKey, h);
    top->fitTTBar    .push_back(hfit[TTBar    ]);
    top->fitHadTop   .push_back(hfit[HadTop   ]);
    top->fitLepTop   .push_back(hfit[LepTop   ]);
    top->fitHadW     .push_back(hfit[HadW     ]);
    top->fitLepW     .push_back(hfit[LepW     ]);
    top->fitHadB     .push_back(hfit[HadB     ]);
    top->fitLightQ   .push_back(hfit[LightQ   ]);
    top->fitLightQBar.push_back(hfit[LightQBar]);
    top->fitLepB     .push_back(hfit[LepB     ]);
    top->fitLepton   .push_back(hfit[Lepton   ]);
    top->fitNeutrino .push_back(hfit[Neutrino ]);

    if(hypoClassKey == TtEvent::kHitFit){
      top->fitProb.push_back(ttEvent->hitFitProb(h));
      top->fitChi2.push_back(ttEvent->hitFitChi2(h));
    }
    else if(hypoClassKey == TtEvent::kKinFit){
      top->fitProb.push_back(ttEvent->fitProb(h));
      top->fitChi2.push_back(ttEvent->fitChi2(h));
    }
    else{
      std::cout << "The given 'hypoClassKey' (" << hypoClassKey << ") is not allowed.\n"
          << "List of allowed 'hypoClassKey':\n"
          << "kinFit: " << TtEvent::kKinFit << "\n"
          << "hitFit: " << TtEvent::kHitFit << std::endl;
    }

    if(semiLepTtEvent){
      std::vector<TLorentzVector> hreco = getPartons(semiLepTtEvent, TtEvent::kMVADisc,
          ttEvent->correspondingHypo(hypoClassKey, h, TtEvent::kMVADisc));
      top->recoTTBar    .push_back(hreco[TTBar    ]);
      top->recoHadTop   .push_back(hreco[HadTop   ]);
      top->recoLepTop   .push_back(hreco[LepTop   ]);
      top->recoHadW     .push_back(hreco[HadW     ]);
      top->recoLepW     .push_back(hreco[LepW     ]);
      top->recoHadB     .push_back(hreco[HadB     ]);
      top->recoLightQ   .push_back(hreco[LightQ   ]);
      top->recoLightQBar.push_back(hreco[LightQBar]);
      top->recoLepB     .push_back(hreco[LepB     ]);
      top->recoLepton   .push_back(hreco[Lepton   ]);
      top->recoNeutrino .push_back(hreco[Neutrino ]);

      top->recoJetIdxHadB.push_back(ttEvent->jetLeptonCombination(hypoClassKey, h)[TtSemiLepEvtPartons::HadB]);
      top->recoJetIdxLightQ.push_back(ttEvent->jetLeptonCombination(hypoClassKey, h)[TtSemiLepEvtPartons::LightQ]);
      top->recoJetIdxLightQBar.push_back(ttEvent->jetLeptonCombination(hypoClassKey, h)[TtSemiLepEvtPartons::LightQBar]);
      top->recoJetIdxLepB.push_back(ttEvent->jetLeptonCombination(hypoClassKey, h)[TtSemiLepEvtPartons::LepB]);
      top->recoJetIdxLepton.push_back(0);
      top->recoJetIdxNeutrino.push_back(0);
    }

    // Fetch reconstructed permutations without fit solution

    for (unsigned int h = 0; h < ttEvent->numberOfAvailableHypos(TtEvent::kMVADisc); ++h) {
      if (!ttEvent->isHypoValid(TtEvent::kMVADisc, h)) break;
      if ( ttEvent->correspondingHypo(TtEvent::kMVADisc, h, hypoClassKey) > -1) continue;

      std::vector<TLorentzVector> hreco = getPartons(semiLepTtEvent, TtEvent::kMVADisc, h);
      top->recoTTBar    .push_back(hreco[TTBar    ]);
      top->recoHadTop   .push_back(hreco[HadTop   ]);
      top->recoLepTop   .push_back(hreco[LepTop   ]);
      top->recoHadW     .push_back(hreco[HadW     ]);
      top->recoLepW     .push_back(hreco[LepW     ]);
      top->recoHadB     .push_back(hreco[HadB     ]);
      top->recoLightQ   .push_back(hreco[LightQ   ]);
      top->recoLightQBar.push_back(hreco[LightQBar]);
      top->recoLepB     .push_back(hreco[LepB     ]);
      top->recoLepton   .push_back(hreco[Lepton   ]);
      top->recoNeutrino .push_back(hreco[Neutrino ]);

      top->recoJetIdxHadB.push_back(ttEvent->jetLeptonCombination(TtEvent::kMVADisc, h)[TtSemiLepEvtPartons::HadB]);
      top->recoJetIdxLightQ.push_back(ttEvent->jetLeptonCombination(TtEvent::kMVADisc, h)[TtSemiLepEvtPartons::LightQ]);
      top->recoJetIdxLightQBar.push_back(ttEvent->jetLeptonCombination(TtEvent::kMVADisc, h)[TtSemiLepEvtPartons::LightQBar]);
      top->recoJetIdxLepB.push_back(ttEvent->jetLeptonCombination(TtEvent::kMVADisc, h)[TtSemiLepEvtPartons::LepB]);
      top->recoJetIdxLepton.push_back(0);
      top->recoJetIdxNeutrino.push_back(0);
    }
  }
  trs->Fill();
}

void 
EventHypothesisAnalyzer::beginJob()
{
  if( !trs ) throw edm::Exception( edm::errors::Configuration, "TreeRegistryService is not registered in cfg file!" );

  top = new TopEvent();
  trs->Branch("top", top);
}

std::vector<TLorentzVector> EventHypothesisAnalyzer::getPartons(const TtSemiLeptonicEvent *ttEvent, TtEvent::HypoClassKey hypoClassKey, unsigned int h) {
  // TODO: validate hypo, catch h=-1
  std::vector<TLorentzVector> parton;

  parton.push_back(TLorentzVector(
      ttEvent->hadronicDecayTop(hypoClassKey, h)->px()     + ttEvent->leptonicDecayTop(hypoClassKey, h)->px(),
      ttEvent->hadronicDecayTop(hypoClassKey, h)->py()     + ttEvent->leptonicDecayTop(hypoClassKey, h)->py(),
      ttEvent->hadronicDecayTop(hypoClassKey, h)->pz()     + ttEvent->leptonicDecayTop(hypoClassKey, h)->pz(),
      ttEvent->hadronicDecayTop(hypoClassKey, h)->energy() + ttEvent->leptonicDecayTop(hypoClassKey, h)->energy()));

  parton.push_back(TLorentzVector(
      ttEvent->hadronicDecayTop(hypoClassKey, h)->px(), ttEvent->hadronicDecayTop(hypoClassKey, h)->py(),
      ttEvent->hadronicDecayTop(hypoClassKey, h)->pz(), ttEvent->hadronicDecayTop(hypoClassKey, h)->energy()));
  parton.push_back(TLorentzVector(
      ttEvent->leptonicDecayTop(hypoClassKey, h)->px(), ttEvent->leptonicDecayTop(hypoClassKey, h)->py(),
      ttEvent->leptonicDecayTop(hypoClassKey, h)->pz(), ttEvent->leptonicDecayTop(hypoClassKey, h)->energy()));

  parton.push_back(TLorentzVector(
      ttEvent->hadronicDecayW(hypoClassKey, h)->px(), ttEvent->hadronicDecayW(hypoClassKey, h)->py(),
      ttEvent->hadronicDecayW(hypoClassKey, h)->pz(), ttEvent->hadronicDecayW(hypoClassKey, h)->energy()));
  parton.push_back(TLorentzVector(
      ttEvent->leptonicDecayW(hypoClassKey, h)->px(), ttEvent->leptonicDecayW(hypoClassKey, h)->py(),
      ttEvent->leptonicDecayW(hypoClassKey, h)->pz(), ttEvent->leptonicDecayW(hypoClassKey, h)->energy()));

  parton.push_back(TLorentzVector(
      ttEvent->hadronicDecayB(hypoClassKey, h)->px(), ttEvent->hadronicDecayB(hypoClassKey, h)->py(),
      ttEvent->hadronicDecayB(hypoClassKey, h)->pz(), ttEvent->hadronicDecayB(hypoClassKey, h)->energy()));
  parton.push_back(TLorentzVector(
      ttEvent->hadronicDecayQuark(hypoClassKey, h)->px(), ttEvent->hadronicDecayQuark(hypoClassKey, h)->py(),
      ttEvent->hadronicDecayQuark(hypoClassKey, h)->pz(), ttEvent->hadronicDecayQuark(hypoClassKey, h)->energy()));
  parton.push_back(TLorentzVector(
      ttEvent->hadronicDecayQuarkBar(hypoClassKey, h)->px(), ttEvent->hadronicDecayQuarkBar(hypoClassKey, h)->py(),
      ttEvent->hadronicDecayQuarkBar(hypoClassKey, h)->pz(), ttEvent->hadronicDecayQuarkBar(hypoClassKey, h)->energy()));
  parton.push_back(TLorentzVector(
      ttEvent->leptonicDecayB(hypoClassKey, h)->px(), ttEvent->leptonicDecayB(hypoClassKey, h)->py(),
      ttEvent->leptonicDecayB(hypoClassKey, h)->pz(), ttEvent->leptonicDecayB(hypoClassKey, h)->energy()));
  parton.push_back(TLorentzVector(
      ttEvent->singleLepton(hypoClassKey, h)->px(),  ttEvent->singleLepton(hypoClassKey, h)->py(),
      ttEvent->singleLepton(hypoClassKey, h)->pz(),  ttEvent->singleLepton(hypoClassKey, h)->energy()));
  parton.push_back(TLorentzVector(
      ttEvent->singleNeutrino(hypoClassKey, h)->px(), ttEvent->singleNeutrino(hypoClassKey, h)->py(),
      ttEvent->singleNeutrino(hypoClassKey, h)->pz(), ttEvent->singleNeutrino(hypoClassKey, h)->energy()));

  return parton;
}

std::vector<TLorentzVector> EventHypothesisAnalyzer::getPartons(const TtFullHadronicEvent *ttEvent, TtEvent::HypoClassKey hypoClassKey, unsigned int h) {
  // TODO: validate hypo, catch h=-1
  std::vector<TLorentzVector> parton;

  parton.push_back(TLorentzVector(
      ttEvent->top(hypoClassKey, h)->px()     + ttEvent->topBar(hypoClassKey, h)->px(),
      ttEvent->top(hypoClassKey, h)->py()     + ttEvent->topBar(hypoClassKey, h)->py(),
      ttEvent->top(hypoClassKey, h)->pz()     + ttEvent->topBar(hypoClassKey, h)->pz(),
      ttEvent->top(hypoClassKey, h)->energy() + ttEvent->topBar(hypoClassKey, h)->energy()));

  parton.push_back(TLorentzVector(
      ttEvent->top(hypoClassKey, h)->px(), ttEvent->top(hypoClassKey, h)->py(),
      ttEvent->top(hypoClassKey, h)->pz(), ttEvent->top(hypoClassKey, h)->energy()));
  parton.push_back(TLorentzVector(
      ttEvent->topBar(hypoClassKey, h)->px(), ttEvent->topBar(hypoClassKey, h)->py(),
      ttEvent->topBar(hypoClassKey, h)->pz(), ttEvent->topBar(hypoClassKey, h)->energy()));

  parton.push_back(TLorentzVector(
      ttEvent->wPlus(hypoClassKey, h)->px(), ttEvent->wPlus(hypoClassKey, h)->py(),
      ttEvent->wPlus(hypoClassKey, h)->pz(), ttEvent->wPlus(hypoClassKey, h)->energy()));
  parton.push_back(TLorentzVector(
      ttEvent->wMinus(hypoClassKey, h)->px(), ttEvent->wMinus(hypoClassKey, h)->py(),
      ttEvent->wMinus(hypoClassKey, h)->pz(), ttEvent->wMinus(hypoClassKey, h)->energy()));

  parton.push_back(TLorentzVector(
      ttEvent->b(hypoClassKey, h)->px(), ttEvent->b(hypoClassKey, h)->py(),
      ttEvent->b(hypoClassKey, h)->pz(), ttEvent->b(hypoClassKey, h)->energy()));
  parton.push_back(TLorentzVector(
      ttEvent->lightQ(hypoClassKey, h)->px(), ttEvent->lightQ(hypoClassKey, h)->py(),
      ttEvent->lightQ(hypoClassKey, h)->pz(), ttEvent->lightQ(hypoClassKey, h)->energy()));
  parton.push_back(TLorentzVector(
      ttEvent->lightQBar(hypoClassKey, h)->px(), ttEvent->lightQBar(hypoClassKey, h)->py(),
      ttEvent->lightQBar(hypoClassKey, h)->pz(), ttEvent->lightQBar(hypoClassKey, h)->energy()));
  parton.push_back(TLorentzVector(
      ttEvent->bBar(hypoClassKey, h)->px(), ttEvent->bBar(hypoClassKey, h)->py(),
      ttEvent->bBar(hypoClassKey, h)->pz(), ttEvent->bBar(hypoClassKey, h)->energy()));
  parton.push_back(TLorentzVector(
      ttEvent->lightP(hypoClassKey, h)->px(), ttEvent->lightP(hypoClassKey, h)->py(),
      ttEvent->lightP(hypoClassKey, h)->pz(), ttEvent->lightP(hypoClassKey, h)->energy()));
  parton.push_back(TLorentzVector(
      ttEvent->lightPBar(hypoClassKey, h)->px(), ttEvent->lightPBar(hypoClassKey, h)->py(),
      ttEvent->lightPBar(hypoClassKey, h)->pz(), ttEvent->lightPBar(hypoClassKey, h)->energy()));

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
