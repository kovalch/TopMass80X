#include "TopMass/TopEventTree/interface/TopEvent.h"
#include "TopMass/TopEventTree/interface/TopEvent_LinkDef.h"

TopEvent::TopEvent() : TObject()
{
  init();
}

void TopEvent::init()
{
  run = -1; lumiBlock = -1; event = -1;
  
  decayChannel = -10;
  
  genpartonTTBar = TLorentzVector();
  genpartonTop1  = TLorentzVector(); genpartonTop2    = TLorentzVector();
  genpartonW1    = TLorentzVector(); genpartonW2      = TLorentzVector();
  genpartonB1    = TLorentzVector(); genpartonW1Prod1 = TLorentzVector(); genpartonW1Prod2 = TLorentzVector();
  genpartonB2    = TLorentzVector(); genpartonW2Prod1 = TLorentzVector(); genpartonW2Prod2 = TLorentzVector();
  
  genpartonJetIdxB1 = -1; genpartonJetIdxW1Prod1 = -1; genpartonJetIdxW1Prod2 = -1;
  genpartonJetIdxB2 = -1; genpartonJetIdxW2Prod1 = -1; genpartonJetIdxW2Prod2 = -1;
  
  recoTTBar  .clear();
  recoTop1   .clear();
  recoTop2   .clear();
  recoW1     .clear();
  recoW2     .clear();
  recoB1     .clear();
  recoW1Prod1.clear();
  recoW1Prod2.clear();
  recoB2     .clear();
  recoW2Prod1.clear();
  recoW2Prod2.clear();

  recoJetIdxB1     .clear();
  recoJetIdxW1Prod1.clear();
  recoJetIdxW1Prod2.clear();
  recoJetIdxB2     .clear();
  recoJetIdxW2Prod1.clear();
  recoJetIdxW2Prod2.clear();
  
  fitTTBar    .clear();
  fitTop1   .clear();
  fitTop2   .clear();
  fitW1     .clear();
  fitW2     .clear();
  fitB1     .clear();
  fitW1Prod1.clear();
  fitW1Prod2.clear();
  fitB2     .clear();
  fitW2Prod1.clear();
  fitW2Prod2.clear();

  combinationType.clear();

  fitProb .clear();
  fitChi2 .clear();
  fitSigMT.clear();
}
