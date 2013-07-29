#include "TopMass/TopEventTree/interface/TopEvent.h"

TopEvent::TopEvent() {
  genpartonTop1       = &genpartonHadTop;
  genpartonTop2       = &genpartonLepTop;
  genpartonW1         = &genpartonHadW;
  genpartonW2         = &genpartonLepW;
  genpartonB1         = &genpartonHadB;
  genpartonLightQ1    = &genpartonLightQ;
  genpartonLightQBar1 = &genpartonLightQBar;
  genpartonB2         = &genpartonLepB;
  genpartonLightQ2    = &genpartonLepton;
  genpartonLightQBar2 = &genpartonNeutrino;
  
  recoTop1       = &recoHadTop;
  recoTop2       = &recoLepTop;
  recoW1         = &recoHadW;
  recoW2         = &recoLepW;
  recoB1         = &recoHadB;
  recoLightQ1    = &recoLightQ;
  recoLightQBar1 = &recoLightQBar;
  recoB2         = &recoLepB;
  recoLightQ2    = &recoLepton;
  recoLightQBar2 = &recoNeutrino;

  recoJetIdxB1         = &recoJetIdxHadB;
  recoJetIdxLightQ1    = &recoJetIdxLightQ;
  recoJetIdxLightQBar1 = &recoJetIdxLightQBar;
  recoJetIdxB2         = &recoJetIdxLepB;
  recoJetIdxLightQ2    = &recoJetIdxLepton;
  recoJetIdxLightQBar2 = &recoJetIdxNeutrino;

  fitTop1       = &fitHadTop;
  fitTop2       = &fitLepTop;
  fitW1         = &fitHadW;
  fitW2         = &fitLepW;
  fitB1         = &fitHadB;
  fitLightQ1    = &fitLightQ;
  fitLightQBar1 = &fitLightQBar;
  fitB2         = &fitLepB;
  fitLightQ2    = &fitLepton;
  fitLightQBar2 = &fitNeutrino;

  init();
}

void TopEvent::init() {
  run = -1; lumiBlock = -1; event = -1;
  
  decayChannel = -10;
  
  genpartonHadTop = TLorentzVector(); genpartonLepTop = TLorentzVector();
  genpartonHadW = TLorentzVector(); genpartonLepW = TLorentzVector(); 
  genpartonHadB = TLorentzVector(); genpartonLightQ = TLorentzVector(); genpartonLightQBar = TLorentzVector();
  genpartonLepB = TLorentzVector(); genpartonLepton = TLorentzVector(); genpartonNeutrino = TLorentzVector();
  
  genpartonJetIdxHadB = -1; genpartonJetIdxLightQ = -1; genpartonJetIdxLightQBar = -1;
  genpartonJetIdxLepB = -1; genpartonJetIdxLepton = -1; genpartonJetIdxNeutrino = -1;
  
  genpartonJetIdxB1 = -1; genpartonJetIdxLightQ1 = -1; genpartonJetIdxLightQBar1 = -1;
  genpartonJetIdxB2 = -1; genpartonJetIdxLightQ2 = -1; genpartonJetIdxLightQBar2 = -1;
  
  recoTTBar    .clear();
  recoHadTop   .clear();
  recoLepTop   .clear();
  recoHadW     .clear();
  recoLepW     .clear();
  recoHadB     .clear();
  recoLightQ   .clear();
  recoLightQBar.clear();
  recoLepB     .clear();
  recoLepton   .clear();
  recoNeutrino .clear();
  recoJetIdxHadB     .clear();
  recoJetIdxLightQ   .clear();
  recoJetIdxLightQBar.clear();
  recoJetIdxLepB     .clear();
  recoJetIdxLepton   .clear();
  recoJetIdxNeutrino .clear();
  
  fitTTBar    .clear();
  fitHadTop   .clear();
  fitLepTop   .clear();
  fitHadW     .clear();
  fitLepW     .clear();
  fitHadB     .clear();
  fitLightQ   .clear();
  fitLightQBar.clear();
  fitLepB     .clear();
  fitLepton   .clear();
  fitNeutrino .clear();

  combinationType.clear();

  fitProb.clear();
  fitChi2.clear();
  fitSigMT.clear();
}

