//#include <memory>

// user include files
//#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TopMass/TopEventTree/interface/TreeRegistryService.h"

#include "AnalysisDataFormats/TopObjects/interface/TtSemiLepEvtPartons.h"
#include "AnalysisDataFormats/TopObjects/interface/TtFullHadEvtPartons.h"
//#include "DataFormats/PatCandidates/interface/Muon.h"

#include "TopMass/TopEventTree/plugins/EventHypothesisAnalyzer.h"

EventHypothesisAnalyzer::EventHypothesisAnalyzer(const edm::ParameterSet& cfg):
ttEvent_     (cfg.getParameter<edm::InputTag>("ttEvent")),
hypoClassKey_(cfg.getParameter<edm::InputTag>("hypoClassKey")),
ttEventGen2_ (cfg.getParameter<edm::InputTag>("ttEventGen2")),

jets_        (cfg.getParameter<edm::InputTag>("jets")), // needed in fullHad channel for reco masses
lepton_      (cfg.getParameter<int>("lepton")),
//leps_        (cfg.getParameter<edm::InputTag>("leps")),
//mets_        (cfg.getParameter<edm::InputTag>("mets")),

//kJetMAX_(cfg.getParameter<int>("maxNJets")),
kMAXCombo_(cfg.getParameter<int>("maxCombo")),
top(0)
{
}

void
EventHypothesisAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& setup)
{
  //////////////////////////////////////////////////////////////////////////////
  // INIT TopEvent
  ////////////////////////////////////////////////////////////////////////////

  top->init();

  top->run       = evt.eventAuxiliary().run()            ;
  top->lumiBlock = evt.eventAuxiliary().luminosityBlock();
  top->event     = evt.eventAuxiliary().event()          ;

  //////////////////////////////////////////////////////////////////////////////
  // get a handle for
  // TtSemiLeptonicEvent / TtFullHadronicEvent
  // and a key to the hypothesis
  //////////////////////////////////////////////////////////////////////////

  edm::Handle<TtSemiLeptonicEvent> hSemiLepTtEvent;
  edm::Handle<TtFullHadronicEvent> hFullHadTtEvent;
  const TtSemiLeptonicEvent *semiLepTtEvent = 0;
  const TtFullHadronicEvent *fullHadTtEvent = 0;
  const TtEvent *ttEvent = 0;

  // needed in fullHad channel for reco masses
  edm::Handle<std::vector<pat::Jet> > hJets;
  const std::vector<pat::Jet> *jets = 0;

  if(ttEvent_.label().find("SemiLep")!=std::string::npos){
    evt.getByLabel(ttEvent_, hSemiLepTtEvent);
    semiLepTtEvent = hSemiLepTtEvent.product();
    ttEvent = semiLepTtEvent;
  }
  else if(ttEvent_.label().find("FullHad")!=std::string::npos){
    evt.getByLabel(ttEvent_, hFullHadTtEvent);
    fullHadTtEvent = hFullHadTtEvent.product();
    ttEvent = fullHadTtEvent;
    evt.getByLabel(jets_, hJets);
    jets = hJets.product();
  }
  else{
    std::cout << "The given 'ttEvent' label (" << ttEvent_.label() << ") is not allowed.\n"
        << "It has to contain either 'SemiLep' or 'FullHad'!" << std::endl;
  }

  top->decayChannel = fillGenPartons(ttEvent->genEvent().get());

  edm::Handle<int> hypoClassKeyHandle;
  evt.getByLabel(hypoClassKey_, hypoClassKeyHandle);

  TtEvent::HypoClassKey& hypoClassKey = (TtEvent::HypoClassKey&) *hypoClassKeyHandle;

  //////////////////////////////////////////////////////////////////////////////
  // Ratteldiekatz (c) Martin Görner
  ////////////////////////////////////////////////////////////////////////////
  // Sorry Markus, mir ist da gerade nichts eingefallen,
  // daher habe ich einen Kommentar dazu von Martin eingefügt
  /////////////////////////////////////////////////////////////////////////

  if(semiLepTtEvent){
    if (semiLepTtEvent->isHypoValid(TtEvent::kGenMatch)) {
      top->genpartonJetIdxB1      = semiLepTtEvent->jetLeptonCombination(TtEvent::kGenMatch)[TtSemiLepEvtPartons::HadB     ];
      top->genpartonJetIdxW1Prod1 = semiLepTtEvent->jetLeptonCombination(TtEvent::kGenMatch)[TtSemiLepEvtPartons::LightQ   ];
      top->genpartonJetIdxW1Prod2 = semiLepTtEvent->jetLeptonCombination(TtEvent::kGenMatch)[TtSemiLepEvtPartons::LightQBar];
      top->genpartonJetIdxB2      = semiLepTtEvent->jetLeptonCombination(TtEvent::kGenMatch)[TtSemiLepEvtPartons::LepB     ];
    }
  }
  else if(fullHadTtEvent){
    if (fullHadTtEvent->isHypoValid(TtEvent::kGenMatch)) {
      top->genpartonJetIdxB1      = fullHadTtEvent->jetLeptonCombination(TtEvent::kGenMatch)[TtFullHadEvtPartons::B        ];
      top->genpartonJetIdxW1Prod1 = fullHadTtEvent->jetLeptonCombination(TtEvent::kGenMatch)[TtFullHadEvtPartons::LightQ   ];
      top->genpartonJetIdxW1Prod2 = fullHadTtEvent->jetLeptonCombination(TtEvent::kGenMatch)[TtFullHadEvtPartons::LightQBar];
      top->genpartonJetIdxB2      = fullHadTtEvent->jetLeptonCombination(TtEvent::kGenMatch)[TtFullHadEvtPartons::BBar     ];
      top->genpartonJetIdxW2Prod1 = fullHadTtEvent->jetLeptonCombination(TtEvent::kGenMatch)[TtFullHadEvtPartons::LightP   ];
      top->genpartonJetIdxW2Prod2 = fullHadTtEvent->jetLeptonCombination(TtEvent::kGenMatch)[TtFullHadEvtPartons::LightPBar];
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // FIT & RECO
  ////////////////////////////////////////////////////////////////////////////

  for (unsigned int h = 0, maxHypos = std::min(ttEvent->numberOfAvailableHypos(hypoClassKey), kMAXCombo_); h < maxHypos; ++h) {
    if (!ttEvent->isHypoValid(hypoClassKey, h)) break;
    std::vector<TLorentzVector> hfit;
    if     (semiLepTtEvent) hfit = getPartons(semiLepTtEvent, hypoClassKey, h);
    else if(fullHadTtEvent) hfit = getPartons(fullHadTtEvent, hypoClassKey, h);
    top->fitTTBar  .push_back(hfit[TopEvent::TTBar    ]);
    top->fitTop1   .push_back(hfit[TopEvent::HadTop   ]);
    top->fitTop2   .push_back(hfit[TopEvent::LepTop   ]);
    top->fitW1     .push_back(hfit[TopEvent::HadW     ]);
    top->fitW2     .push_back(hfit[TopEvent::LepW     ]);
    top->fitB1     .push_back(hfit[TopEvent::HadB     ]);
    top->fitW1Prod1.push_back(hfit[TopEvent::LightQ   ]);
    top->fitW1Prod2.push_back(hfit[TopEvent::LightQBar]);
    top->fitB2     .push_back(hfit[TopEvent::LepB     ]);
    top->fitW2Prod1.push_back(hfit[TopEvent::Lepton   ]);
    top->fitW2Prod2.push_back(hfit[TopEvent::Neutrino ]);

    if(hypoClassKey == TtEvent::kHitFit){
      top->fitProb .push_back(ttEvent->hitFitProb(h));
      top->fitChi2 .push_back(ttEvent->hitFitChi2(h));
      top->fitSigMT.push_back(ttEvent->hitFitSigMT(h));
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
      top->recoJetIdxB1     .push_back(ttEvent->jetLeptonCombination(hypoClassKey, h)[TtSemiLepEvtPartons::HadB     ]);
      top->recoJetIdxW1Prod1.push_back(ttEvent->jetLeptonCombination(hypoClassKey, h)[TtSemiLepEvtPartons::LightQ   ]);
      top->recoJetIdxW1Prod2.push_back(ttEvent->jetLeptonCombination(hypoClassKey, h)[TtSemiLepEvtPartons::LightQBar]);
      top->recoJetIdxB2     .push_back(ttEvent->jetLeptonCombination(hypoClassKey, h)[TtSemiLepEvtPartons::LepB     ]);
      top->recoJetIdxW2Prod1.push_back(0);
      top->recoJetIdxW2Prod2.push_back(0);

      std::vector<TLorentzVector> hreco = getPartons(semiLepTtEvent, TtEvent::kMVADisc,
          ttEvent->correspondingHypo(hypoClassKey, h, TtEvent::kMVADisc));
      top->recoTTBar  .push_back(hreco[TopEvent::TTBar    ]);
      top->recoTop1   .push_back(hreco[TopEvent::HadTop   ]);
      top->recoTop2   .push_back(hreco[TopEvent::LepTop   ]);
      top->recoW1     .push_back(hreco[TopEvent::HadW     ]);
      top->recoW2     .push_back(hreco[TopEvent::LepW     ]);
      top->recoB1     .push_back(hreco[TopEvent::HadB     ]);
      top->recoW1Prod1.push_back(hreco[TopEvent::LightQ   ]);
      top->recoW1Prod2.push_back(hreco[TopEvent::LightQBar]);
      top->recoB2     .push_back(hreco[TopEvent::LepB     ]);
      top->recoW2Prod1.push_back(hreco[TopEvent::Lepton   ]);
      top->recoW2Prod2.push_back(hreco[TopEvent::Neutrino ]);
      
      if(semiLepTtEvent->isHypoValid(TtEvent::kGenMatch)){
        top->combinationType.push_back(comboTypeSemiLep());
      }
      else{
        top->combinationType.push_back(6);
      }
    }
    else if(fullHadTtEvent){
      top->recoJetIdxB1     .push_back(ttEvent->jetLeptonCombination(hypoClassKey, h)[TtFullHadEvtPartons::B        ]);
      top->recoJetIdxW1Prod1.push_back(ttEvent->jetLeptonCombination(hypoClassKey, h)[TtFullHadEvtPartons::LightQ   ]);
      top->recoJetIdxW1Prod2.push_back(ttEvent->jetLeptonCombination(hypoClassKey, h)[TtFullHadEvtPartons::LightQBar]);
      top->recoJetIdxB2     .push_back(ttEvent->jetLeptonCombination(hypoClassKey, h)[TtFullHadEvtPartons::BBar     ]);
      top->recoJetIdxW2Prod1.push_back(ttEvent->jetLeptonCombination(hypoClassKey, h)[TtFullHadEvtPartons::LightP   ]);
      top->recoJetIdxW2Prod2.push_back(ttEvent->jetLeptonCombination(hypoClassKey, h)[TtFullHadEvtPartons::LightPBar]);

      // do this AFTER *recoJetIdx*'s have been filled, values are used here !!!
      std::vector<TLorentzVector> hreco = getPartons(jets);
      top->recoTTBar  .push_back(hreco[TopEvent::TTBar    ]);
      top->recoTop1   .push_back(hreco[TopEvent::HadTop   ]);
      top->recoTop2   .push_back(hreco[TopEvent::LepTop   ]);
      top->recoW1     .push_back(hreco[TopEvent::HadW     ]);
      top->recoW2     .push_back(hreco[TopEvent::LepW     ]);
      top->recoB1     .push_back(hreco[TopEvent::HadB     ]);
      top->recoW1Prod1.push_back(hreco[TopEvent::LightQ   ]);
      top->recoW1Prod2.push_back(hreco[TopEvent::LightQBar]);
      top->recoB2     .push_back(hreco[TopEvent::LepB     ]);
      top->recoW2Prod1.push_back(hreco[TopEvent::Lepton   ]);
      top->recoW2Prod2.push_back(hreco[TopEvent::Neutrino ]);

      // get 2. genMatch to see if the case is ambiguous or unmatchable
      edm::Handle<TtFullHadronicEvent> hFullHadTtEvent2;
      evt.getByLabel(ttEventGen2_, hFullHadTtEvent2);
      bool genMatch1Valid = fullHadTtEvent->isHypoValid(TtEvent::kGenMatch);
      bool genMatch2Valid = (hFullHadTtEvent2.isValid() && hFullHadTtEvent2->isHypoValid(TtEvent::kGenMatch));

      if( genMatch1Valid ){
        // do this AFTER *recoJetIdx* have been filled, values are used here !!!
        top->combinationType.push_back(comboTypeFullHad());
      }
      else{
        if( genMatch2Valid ){
          top->combinationType.push_back(5);
        }
        else{
          top->combinationType.push_back(6);
        }
      }
    }
    if(semiLepTtEvent){
      //////////////////////////////////////////////////////////////////////////////
      // Fetch reconstructed permutations without fit solution, get leptonFlavour
      // not implemented for fullHad version up to now
      ///////////////////////////////////////////////////////////////////////////
      
      top->leptonFlavour = -lepton_ * semiLepTtEvent->singleLepton(TtEvent::kMVADisc)->charge();
      
      std::cout << "leptonFlavour: " << top->leptonFlavour << std::endl;
      
      for (unsigned int h = 0; h < ttEvent->numberOfAvailableHypos(TtEvent::kMVADisc); ++h) {
        if (!ttEvent->isHypoValid(TtEvent::kMVADisc, h)) break;
        if ( ttEvent->correspondingHypo(TtEvent::kMVADisc, h, hypoClassKey) > -1) continue;

        std::vector<TLorentzVector> hreco = getPartons(semiLepTtEvent, TtEvent::kMVADisc, h);
        top->recoTTBar  .push_back(hreco[TopEvent::TTBar    ]);
        top->recoTop1   .push_back(hreco[TopEvent::HadTop   ]);
        top->recoTop2   .push_back(hreco[TopEvent::LepTop   ]);
        top->recoW1     .push_back(hreco[TopEvent::HadW     ]);
        top->recoW2     .push_back(hreco[TopEvent::LepW     ]);
        top->recoB1     .push_back(hreco[TopEvent::HadB     ]);
        top->recoW1Prod1.push_back(hreco[TopEvent::LightQ   ]);
        top->recoW1Prod2.push_back(hreco[TopEvent::LightQBar]);
        top->recoB2     .push_back(hreco[TopEvent::LepB     ]);
        top->recoW2Prod1.push_back(hreco[TopEvent::Lepton   ]);
        top->recoW2Prod2.push_back(hreco[TopEvent::Neutrino ]);

        top->recoJetIdxB1     .push_back(ttEvent->jetLeptonCombination(TtEvent::kMVADisc, h)[TtSemiLepEvtPartons::HadB]);
        top->recoJetIdxW1Prod1.push_back(ttEvent->jetLeptonCombination(TtEvent::kMVADisc, h)[TtSemiLepEvtPartons::LightQ]);
        top->recoJetIdxW1Prod2.push_back(ttEvent->jetLeptonCombination(TtEvent::kMVADisc, h)[TtSemiLepEvtPartons::LightQBar]);
        top->recoJetIdxB2     .push_back(ttEvent->jetLeptonCombination(TtEvent::kMVADisc, h)[TtSemiLepEvtPartons::LepB]);
        top->recoJetIdxW2Prod1.push_back(0);
        top->recoJetIdxW2Prod2.push_back(0);
      }
    }
  }
  trs->Fill();
}

void 
EventHypothesisAnalyzer::beginJob()
{
  if( !trs ) throw edm::Exception( edm::errors::Configuration, "TreeRegistryService is not registered in cfg file!" );

  top = new TopEvent();
  trs->Branch("top.", top);
}

std::vector<TLorentzVector> EventHypothesisAnalyzer::getPartons(const TtSemiLeptonicEvent *ttEvent, TtEvent::HypoClassKey hypoClassKey, unsigned int h)
{
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

std::vector<TLorentzVector> EventHypothesisAnalyzer::getPartons(const TtFullHadronicEvent *ttEvent, TtEvent::HypoClassKey hypoClassKey, unsigned int h)
{
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

std::vector<TLorentzVector> EventHypothesisAnalyzer::getPartons(const std::vector<pat::Jet> *jets)
{
  int idxB1         = top->recoJetIdxB1     .back();
  int idxLightQ1    = top->recoJetIdxW1Prod1.back();
  int idxLightQBar1 = top->recoJetIdxW1Prod2.back();
  int idxB2         = top->recoJetIdxB2     .back();
  int idxLightQ2    = top->recoJetIdxW2Prod1.back();
  int idxLightQBar2 = top->recoJetIdxW2Prod2.back();

  // TODO: validate hypo, catch h=-1
  std::vector<TLorentzVector> parton;

  parton.push_back(TLorentzVector(
      jets->at(idxB1).px()     + jets->at(idxLightQ1).px()     + jets->at(idxLightQBar1).px()     +  jets->at(idxB2).px()     + jets->at(idxLightQ2).px()     + jets->at(idxLightQBar2).px(),
      jets->at(idxB1).py()     + jets->at(idxLightQ1).py()     + jets->at(idxLightQBar1).py()     +  jets->at(idxB2).py()     + jets->at(idxLightQ2).py()     + jets->at(idxLightQBar2).py(),
      jets->at(idxB1).pz()     + jets->at(idxLightQ1).pz()     + jets->at(idxLightQBar1).pz()     +  jets->at(idxB2).pz()     + jets->at(idxLightQ2).pz()     + jets->at(idxLightQBar2).pz(),
      jets->at(idxB1).energy() + jets->at(idxLightQ1).energy() + jets->at(idxLightQBar1).energy() +  jets->at(idxB2).energy() + jets->at(idxLightQ2).energy() + jets->at(idxLightQBar2).energy()));

  parton.push_back(TLorentzVector(
      jets->at(idxB1).px()     + jets->at(idxLightQ1).px()     + jets->at(idxLightQBar1).px(),
      jets->at(idxB1).py()     + jets->at(idxLightQ1).py()     + jets->at(idxLightQBar1).py(),
      jets->at(idxB1).pz()     + jets->at(idxLightQ1).pz()     + jets->at(idxLightQBar1).pz(),
      jets->at(idxB1).energy() + jets->at(idxLightQ1).energy() + jets->at(idxLightQBar1).energy()));
  parton.push_back(TLorentzVector(
      jets->at(idxB2).px()     + jets->at(idxLightQ2).px()     + jets->at(idxLightQBar2).px(),
      jets->at(idxB2).py()     + jets->at(idxLightQ2).py()     + jets->at(idxLightQBar2).py(),
      jets->at(idxB2).pz()     + jets->at(idxLightQ2).pz()     + jets->at(idxLightQBar2).pz(),
      jets->at(idxB2).energy() + jets->at(idxLightQ2).energy() + jets->at(idxLightQBar2).energy()));

  parton.push_back(TLorentzVector(
      jets->at(idxLightQ1).px() + jets->at(idxLightQBar1).px(), jets->at(idxLightQ1).py()     + jets->at(idxLightQBar1).py(),
      jets->at(idxLightQ1).pz() + jets->at(idxLightQBar1).pz(), jets->at(idxLightQ1).energy() + jets->at(idxLightQBar1).energy()));
  parton.push_back(TLorentzVector(
      jets->at(idxLightQ2).px() + jets->at(idxLightQBar2).px(), jets->at(idxLightQ2).py()     + jets->at(idxLightQBar2).py(),
      jets->at(idxLightQ2).pz() + jets->at(idxLightQBar2).pz(), jets->at(idxLightQ2).energy() + jets->at(idxLightQBar2).energy()));

  parton.push_back(TLorentzVector(jets->at(idxB1        ).px(), jets->at(idxB1        ).py(), jets->at(idxB1        ).pz(), jets->at(idxB1        ).energy()));
  parton.push_back(TLorentzVector(jets->at(idxLightQ1   ).px(), jets->at(idxLightQ1   ).py(), jets->at(idxLightQ1   ).pz(), jets->at(idxLightQ1   ).energy()));
  parton.push_back(TLorentzVector(jets->at(idxLightQBar1).px(), jets->at(idxLightQBar1).py(), jets->at(idxLightQBar1).pz(), jets->at(idxLightQBar1).energy()));
  parton.push_back(TLorentzVector(jets->at(idxB2        ).px(), jets->at(idxB2        ).py(), jets->at(idxB2        ).pz(), jets->at(idxB2        ).energy()));
  parton.push_back(TLorentzVector(jets->at(idxLightQ2   ).px(), jets->at(idxLightQ2   ).py(), jets->at(idxLightQ2   ).pz(), jets->at(idxLightQ2   ).energy()));
  parton.push_back(TLorentzVector(jets->at(idxLightQBar2).px(), jets->at(idxLightQBar2).py(), jets->at(idxLightQBar2).pz(), jets->at(idxLightQBar2).energy()));

  return parton;
}


//////////////////////////////////////////////////////////////////////////////////
// GENPARTON SEMI-LEP
// TTBar, HadTop, LepTop, HadW, LepW,
// HadB, LightQ, LightQBar, LepB, Lepton, Neutrino
//////////////////////////////////////////////////////////////////////////////
// GENPARTON FULL-HAD
// TTBar, Top, TopBar, W+, W-,
// B, LightQ1, LightQBar1, BBar, LightQ2, LightQBar2
//////////////////////////////////////////////////////////////////////////
// GENPARTON FULL-LEP
// TTBar, Top, TopBar, W+, W-,
// B, LeptonBar, Neutrino, BBar, Lepton, NeutrinoBar
//////////////////////////////////////////////////////////////////////
// identify event by the following number coding:
// -10 : undefined
//  -1 : non ttbar MC
//   0 : DATA
//   1 : fully hadronic        ttbar MC
//   2 : semileptonic unknown  ttbar MC
//  21 : semileptonic electron ttbar MC
//  22 : semileptonic muon     ttbar MC
//  23 : semileptonic tau      ttbar MC
//   3 : dileptonic unknown    ttbar MC
// 311 : dileptonic e   e      ttbar MC
// 312 : dileptonic e   mu     ttbar MC
// 313 : dileptonic e   tau    ttbar MC
// 322 : dileptonic mu  mu     ttbar MC
// 323 : dileptonic mu  tau    ttbar MC
// 333 : dileptonic tau tau    ttbar MC
/////////////////////////////////////////////////////

/*

genpartonTTBar;
genpartonHadTop;
genpartonLepTop;
genpartonHadW;
genpartonLepW;
genpartonHadB;
genpartonLightQ;
genpartonLightQBar;
genpartonLepB;
genpartonLepton;
genpartonNeutrino;

*/
int
EventHypothesisAnalyzer::fillGenPartons(const TtGenEvent *genEvent)
{
  int decayChannel = -10;
  // if there is no genEvent classify event as data
  if(      !genEvent            ) decayChannel =  0;
  else if( !genEvent->isTtBar() ) decayChannel = -1;
  // now we are sure to have a ttbar event
  else{
    if(  genEvent->isFullHadronic() ) {
      top->genpartonTTBar.SetPxPyPzE(
          genEvent->top()->px()     + genEvent->topBar()->px(),
          genEvent->top()->py()     + genEvent->topBar()->py(),
          genEvent->top()->pz()     + genEvent->topBar()->pz(),
          genEvent->top()->energy() + genEvent->topBar()->energy());

      top->genpartonTop1.SetPxPyPzE(
          genEvent->top()->px(), genEvent->top()->py(),
          genEvent->top()->pz(), genEvent->top()->energy());
      top->genpartonTop2.SetPxPyPzE(
          genEvent->topBar()->px(), genEvent->topBar()->py(),
          genEvent->topBar()->pz(), genEvent->topBar()->energy());

      top->genpartonW1.SetPxPyPzE(
          genEvent->wPlus()->px(), genEvent->wPlus()->py(),
          genEvent->wPlus()->pz(), genEvent->wPlus()->energy());
      top->genpartonW2.SetPxPyPzE(
          genEvent->wMinus()->px(), genEvent->wMinus()->py(),
          genEvent->wMinus()->pz(), genEvent->wMinus()->energy());

      top->genpartonB1.SetPxPyPzE(
          genEvent->b()->px(), genEvent->b()->py(),
          genEvent->b()->pz(), genEvent->b()->energy());

      top->genpartonW1Prod1.SetPxPyPzE(
          genEvent->daughterQuarkOfWPlus()->px(), genEvent->daughterQuarkOfWPlus()->py(),
          genEvent->daughterQuarkOfWPlus()->pz(), genEvent->daughterQuarkOfWPlus()->energy());
      top->genpartonW1Prod2.SetPxPyPzE(
          genEvent->daughterQuarkBarOfWPlus()->px(), genEvent->daughterQuarkBarOfWPlus()->py(),
          genEvent->daughterQuarkBarOfWPlus()->pz(), genEvent->daughterQuarkBarOfWPlus()->energy());

      top->genpartonB2.SetPxPyPzE(
          genEvent->bBar()->px(), genEvent->bBar()->py(),
          genEvent->bBar()->pz(), genEvent->bBar()->energy());

      top->genpartonW2Prod1.SetPxPyPzE(
          genEvent->daughterQuarkOfWMinus()->px(), genEvent->daughterQuarkOfWMinus()->py(),
          genEvent->daughterQuarkOfWMinus()->pz(), genEvent->daughterQuarkOfWMinus()->energy());
      top->genpartonW2Prod2.SetPxPyPzE(
          genEvent->daughterQuarkBarOfWMinus()->px(), genEvent->daughterQuarkBarOfWMinus()->py(),
          genEvent->daughterQuarkBarOfWMinus()->pz(), genEvent->daughterQuarkBarOfWMinus()->energy());

      decayChannel =  1;
    }
    else if(  genEvent->isSemiLeptonic() ) {
      top->genpartonTTBar.SetPxPyPzE(
          genEvent->hadronicDecayTop()->px()     + genEvent->leptonicDecayTop()->px(),
          genEvent->hadronicDecayTop()->py()     + genEvent->leptonicDecayTop()->py(),
          genEvent->hadronicDecayTop()->pz()     + genEvent->leptonicDecayTop()->pz(),
          genEvent->hadronicDecayTop()->energy() + genEvent->leptonicDecayTop()->energy());

      top->genpartonTop1.SetPxPyPzE(
          genEvent->hadronicDecayTop()->px(), genEvent->hadronicDecayTop()->py(),
          genEvent->hadronicDecayTop()->pz(), genEvent->hadronicDecayTop()->energy());
      top->genpartonTop2.SetPxPyPzE(
          genEvent->leptonicDecayTop()->px(), genEvent->leptonicDecayTop()->py(),
          genEvent->leptonicDecayTop()->pz(), genEvent->leptonicDecayTop()->energy());

      top->genpartonW1.SetPxPyPzE(
          genEvent->hadronicDecayW()->px(), genEvent->hadronicDecayW()->py(),
          genEvent->hadronicDecayW()->pz(), genEvent->hadronicDecayW()->energy());
      top->genpartonW2.SetPxPyPzE(
          genEvent->leptonicDecayW()->px(), genEvent->leptonicDecayW()->py(),
          genEvent->leptonicDecayW()->pz(), genEvent->leptonicDecayW()->energy());

      top->genpartonB1.SetPxPyPzE(
          genEvent->hadronicDecayB()->px(), genEvent->hadronicDecayB()->py(),
          genEvent->hadronicDecayB()->pz(), genEvent->hadronicDecayB()->energy());

      top->genpartonW1Prod1.SetPxPyPzE(
          genEvent->hadronicDecayQuark()->px(), genEvent->hadronicDecayQuark()->py(),
          genEvent->hadronicDecayQuark()->pz(), genEvent->hadronicDecayQuark()->energy());
      top->genpartonW1Prod2.SetPxPyPzE(
          genEvent->hadronicDecayQuarkBar()->px(), genEvent->hadronicDecayQuarkBar()->py(),
          genEvent->hadronicDecayQuarkBar()->pz(), genEvent->hadronicDecayQuarkBar()->energy());

      top->genpartonB2.SetPxPyPzE(
          genEvent->leptonicDecayB()->px(), genEvent->leptonicDecayB()->py(),
          genEvent->leptonicDecayB()->pz(), genEvent->leptonicDecayB()->energy());

      top->genpartonW2Prod1.SetPxPyPzE(
          genEvent->singleLepton()->px(), genEvent->singleLepton()->py(),
          genEvent->singleLepton()->pz(), genEvent->singleLepton()->energy());
      top->genpartonW2Prod2.SetPxPyPzE(
          genEvent->singleNeutrino()->px(), genEvent->singleNeutrino()->py(),
          genEvent->singleNeutrino()->pz(), genEvent->singleNeutrino()->energy());

      switch( genEvent->semiLeptonicChannel() ) {
      case WDecay::kElec : decayChannel = 21; break;
      case WDecay::kMuon : decayChannel = 22; break;
      case WDecay::kTau  : decayChannel = 23; break;
      default            : decayChannel =  2; break;
      }
    }
    else if( genEvent->isFullLeptonic() ) {
      top->genpartonTTBar.SetPxPyPzE(
          genEvent->top()->px()     + genEvent->topBar()->px(),
          genEvent->top()->py()     + genEvent->topBar()->py(),
          genEvent->top()->pz()     + genEvent->topBar()->pz(),
          genEvent->top()->energy() + genEvent->topBar()->energy());

      top->genpartonTop1.SetPxPyPzE(
          genEvent->top()->px(), genEvent->top()->py(),
          genEvent->top()->pz(), genEvent->top()->energy());
      top->genpartonTop2.SetPxPyPzE(
          genEvent->topBar()->px(), genEvent->topBar()->py(),
          genEvent->topBar()->pz(), genEvent->topBar()->energy());

      top->genpartonW1.SetPxPyPzE(
          genEvent->wPlus()->px(), genEvent->wPlus()->py(),
          genEvent->wPlus()->pz(), genEvent->wPlus()->energy());
      top->genpartonW2.SetPxPyPzE(
          genEvent->wMinus()->px(), genEvent->wMinus()->py(),
          genEvent->wMinus()->pz(), genEvent->wMinus()->energy());

      top->genpartonB1.SetPxPyPzE(
          genEvent->b()->px(), genEvent->b()->py(),
          genEvent->b()->pz(), genEvent->b()->energy());

      top->genpartonW1Prod1.SetPxPyPzE(
          genEvent->leptonBar()->px(), genEvent->leptonBar()->py(),
          genEvent->leptonBar()->pz(), genEvent->leptonBar()->energy());
      top->genpartonW1Prod2.SetPxPyPzE(
          genEvent->neutrino()->px(), genEvent->neutrino()->py(),
          genEvent->neutrino()->pz(), genEvent->neutrino()->energy());

      top->genpartonB2.SetPxPyPzE(
          genEvent->bBar()->px(), genEvent->bBar()->py(),
          genEvent->bBar()->pz(), genEvent->bBar()->energy());

      top->genpartonW2Prod1.SetPxPyPzE(
          genEvent->lepton()->px(), genEvent->lepton()->py(),
          genEvent->lepton()->pz(), genEvent->lepton()->energy());
      top->genpartonW2Prod2.SetPxPyPzE(
          genEvent->neutrinoBar()->px(), genEvent->neutrinoBar()->py(),
          genEvent->neutrinoBar()->pz(), genEvent->neutrinoBar()->energy());

      if     (  genEvent->fullLeptonicChannel().first  == WDecay::kElec  &&
                genEvent->fullLeptonicChannel().second == WDecay::kElec  ) decayChannel = 311;
      else if( (genEvent->fullLeptonicChannel().first  == WDecay::kElec  &&
                genEvent->fullLeptonicChannel().second == WDecay::kMuon) ||
               (genEvent->fullLeptonicChannel().first  == WDecay::kMuon  &&
                genEvent->fullLeptonicChannel().second == WDecay::kElec) ) decayChannel = 312;
      else if( (genEvent->fullLeptonicChannel().first  == WDecay::kElec  &&
                genEvent->fullLeptonicChannel().second == WDecay::kTau ) ||
               (genEvent->fullLeptonicChannel().first  == WDecay::kTau   &&
                genEvent->fullLeptonicChannel().second == WDecay::kElec) ) decayChannel = 313;
      else if(  genEvent->fullLeptonicChannel().first  == WDecay::kMuon  &&
                genEvent->fullLeptonicChannel().second == WDecay::kMuon  ) decayChannel = 322;
      else if( (genEvent->fullLeptonicChannel().first  == WDecay::kMuon  &&
                genEvent->fullLeptonicChannel().second == WDecay::kTau ) ||
               (genEvent->fullLeptonicChannel().first  == WDecay::kTau   &&
                genEvent->fullLeptonicChannel().second == WDecay::kMuon) ) decayChannel = 323;
      else if(  genEvent->fullLeptonicChannel().first  == WDecay::kTau   &&
                genEvent->fullLeptonicChannel().second == WDecay::kTau   ) decayChannel = 333;
      else decayChannel = 3;
    }
  }
  return decayChannel;
}

/// function to find types of jet-combinations in KinFits (1 right, 2 one branche right, other branch inner-branch mixup, 3 both branches inner-branch mixup, 4 cross-branch mixup, -1 to -6 number of falsely picked jets)
short
EventHypothesisAnalyzer::comboTypeFullHad()
{
  short comboTypeID = comboTypeIDCalculatorFullHad();
  //std::cout << "ID: " << comboTypeID << std::endl;
  // falsely picked jets
  if(comboTypeID < 0)
    return comboTypeID;

  // correct permutations
  if(comboTypeID ==   0 || comboTypeID ==   1 || comboTypeID ==  24 || comboTypeID ==  25 ||
     comboTypeID == 450 || comboTypeID == 451 || comboTypeID == 474 || comboTypeID == 475)
    return 1;

  // one branch mixup
  if((comboTypeID >=   2 && comboTypeID <=   5) || (comboTypeID >=  26 && comboTypeID <=  29) ||
      comboTypeID == 120 || comboTypeID == 121  ||  comboTypeID == 144 || comboTypeID == 145  ||
      comboTypeID == 240 || comboTypeID == 241  ||  comboTypeID == 264 || comboTypeID == 265  ||
     (comboTypeID >= 452 && comboTypeID <= 455) || (comboTypeID >= 476 && comboTypeID <= 479) ||
      comboTypeID == 570 || comboTypeID == 571  ||  comboTypeID == 594 || comboTypeID == 595  ||
      comboTypeID == 690 || comboTypeID == 691  ||  comboTypeID == 714 || comboTypeID == 715)
    return 2;

  // both branches mixup
  if((comboTypeID >= 122 && comboTypeID <= 125) || (comboTypeID >= 146 && comboTypeID <= 149) ||
     (comboTypeID >= 242 && comboTypeID <= 245) || (comboTypeID >= 266 && comboTypeID <= 269) ||
     (comboTypeID >= 572 && comboTypeID <= 575) || (comboTypeID >= 596 && comboTypeID <= 599) ||
     (comboTypeID >= 692 && comboTypeID <= 695) || (comboTypeID >= 716 && comboTypeID <= 719))
    return 3;

  // mixup cross branches
  if((comboTypeID >=   6 && comboTypeID <=  23) || (comboTypeID >=  30 && comboTypeID <= 119) ||
     (comboTypeID >= 126 && comboTypeID <= 143) || (comboTypeID >= 150 && comboTypeID <= 239) ||
     (comboTypeID >= 246 && comboTypeID <= 263) || (comboTypeID >= 270 && comboTypeID <= 449) ||
     (comboTypeID >= 456 && comboTypeID <= 473) || (comboTypeID >= 480 && comboTypeID <= 569) ||
     (comboTypeID >= 576 && comboTypeID <= 593) || (comboTypeID >= 600 && comboTypeID <= 689) ||
     (comboTypeID >= 696 && comboTypeID <= 713) )
    return 4;

  return -37;
}

short
EventHypothesisAnalyzer::comboTypeSemiLep()
{
  short comboTypeID = comboTypeIDCalculatorSemiLep();
  //std::cout << "ID: " << comboTypeID << std::endl;
  // falsely picked jets
  if(comboTypeID < 0)
    return comboTypeID;

  // correct permutations
  if(comboTypeID ==  0 || comboTypeID ==  2)
    return 1;
  
  // hadronic branch mixup (correct tops)
  if(comboTypeID ==  6 || comboTypeID ==  8 || comboTypeID == 12 || comboTypeID == 14)
    return 2;
  
  // swapped b-quarks (correct W)
  if(comboTypeID == 21 || comboTypeID == 23)
    return 3;
  
  // other mixups (complete mess)
  if(comboTypeID >=  0 && comboTypeID <= 23)
    return 4;
  
  // error
  return -37;
}

/// assign unique ID to every permutation type
short
EventHypothesisAnalyzer::comboTypeIDCalculatorFullHad()
{
  /// vector to store the jet indices
  std::vector<int> jetIndexFit;
  std::vector<int> jetIndexGen;

  jetIndexFit.push_back(top->recoJetIdxB1     .back());
  jetIndexFit.push_back(top->recoJetIdxW1Prod1.back());
  jetIndexFit.push_back(top->recoJetIdxW1Prod2.back());
  jetIndexFit.push_back(top->recoJetIdxB2     .back());
  jetIndexFit.push_back(top->recoJetIdxW2Prod1.back());
  jetIndexFit.push_back(top->recoJetIdxW2Prod2.back());
  
  jetIndexGen.push_back(top->genpartonJetIdxB1     );
  jetIndexGen.push_back(top->genpartonJetIdxW1Prod1);
  jetIndexGen.push_back(top->genpartonJetIdxW1Prod2);
  jetIndexGen.push_back(top->genpartonJetIdxB2     );
  jetIndexGen.push_back(top->genpartonJetIdxW2Prod1);
  jetIndexGen.push_back(top->genpartonJetIdxW2Prod2);
  
  //std::cout << "fit: " << jetIndexFit[0] << " " << jetIndexFit[1] << " " << jetIndexFit[2] << " " << jetIndexFit[3] << " " << jetIndexFit[4] << " " << jetIndexFit[5] << std::endl;
  //std::cout << "gen: " << jetIndexGen[0] << " " << jetIndexGen[1] << " " << jetIndexGen[2] << " " << jetIndexGen[3] << " " << jetIndexGen[4] << " " << jetIndexGen[5] << std::endl;

  return comboTypeAlgo(jetIndexFit, jetIndexGen);
}

short
EventHypothesisAnalyzer::comboTypeIDCalculatorSemiLep()
{
  /// vector to store the jet indices
  std::vector<int> jetIndexFit;
  std::vector<int> jetIndexGen;

  jetIndexFit.push_back(top->recoJetIdxB1     .back());
  jetIndexFit.push_back(top->recoJetIdxW1Prod1.back());
  jetIndexFit.push_back(top->recoJetIdxW1Prod2.back());
  jetIndexFit.push_back(top->recoJetIdxB2     .back());
  
  jetIndexGen.push_back(top->genpartonJetIdxB1     );
  jetIndexGen.push_back(top->genpartonJetIdxW1Prod1);
  jetIndexGen.push_back(top->genpartonJetIdxW1Prod2);
  jetIndexGen.push_back(top->genpartonJetIdxB2     );
  
  return comboTypeAlgo(jetIndexFit, jetIndexGen);
}

short
EventHypothesisAnalyzer::comboTypeAlgo(std::vector<int> jetIndexFit, std::vector<int> jetIndexGen)
{
  short result = 0;
  short wrongJets = 0;
  for(unsigned short iFit = 0; iFit < jetIndexFit.size(); ++iFit){
    short fact = TMath::Factorial(jetIndexFit.size()-iFit-1);
    //std::cout << "result = " << result;
    bool foundPair = false;
    for(unsigned short jGen = 0, jGenCount = 0; jGen < jetIndexGen.size(); ++jGen, ++jGenCount){
      if(jetIndexGen.at(jGen) == -37){
    --jGenCount;
    continue;
      }
      if(jetIndexFit.at(iFit) == jetIndexGen.at(jGen)){
    foundPair = true;
    result += fact*(jGenCount);
    //std::cout << " (" << iFit << "," << jGen << "," << jGenCount << ")";
    jetIndexGen[jGen] = -37;
      }
    }
    if(!foundPair){
      ++wrongJets;
    }
    //std::cout << " -> " << result << std::endl;
  }
  if(wrongJets){
    result = 0-wrongJets;
  }

  //std::cout << "result: " << result << std::endl;

  return result;
}

std::vector<short>
EventHypothesisAnalyzer::comboTypeAlgoInverted(std::vector<int> jetIndexGen, short comboType)
{
  unsigned short lGen = jetIndexGen.size();
  if(comboType<0 || comboType>=TMath::Factorial(lGen)) return std::vector<short>(0);
  std::vector<short> result(lGen);
  for(unsigned short iGen = 0; iGen < lGen; ++iGen){
    short idx = lGen-iGen-1;
    short fact = TMath::Factorial(idx);
    short calc = comboType/fact;
    result[iGen] = jetIndexGen.at(calc);
    jetIndexGen.erase(jetIndexGen.begin()+calc);
    comboType -= calc*fact;
  }
  return result;
}

void
EventHypothesisAnalyzer::endJob() 
{
}

EventHypothesisAnalyzer::~EventHypothesisAnalyzer()
{
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(EventHypothesisAnalyzer);
