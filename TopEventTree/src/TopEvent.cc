#include "TopMass/TopEventTree/interface/TopEvent.h"
#include "TopMass/TopEventTree/interface/TopEvent_LinkDef.h"

TopEvent::TopEvent() : TObject()
{
  init();
}

void TopEvent::init()
{
  _run = -1; _lumiBlock = -1; _event = -1;
  
  _decayChannel = -10;
  
  _genpartonTTBar = TLorentzVector();
  _genpartonTop1  = TLorentzVector(); _genpartonTop2    = TLorentzVector();
  _genpartonW1    = TLorentzVector(); _genpartonW2      = TLorentzVector();
  _genpartonB1    = TLorentzVector(); _genpartonW1Prod1 = TLorentzVector(); _genpartonW1Prod2 = TLorentzVector();
  _genpartonB2    = TLorentzVector(); _genpartonW2Prod1 = TLorentzVector(); _genpartonW2Prod2 = TLorentzVector();
  
  _genpartonJetIdxB1 = -1; _genpartonJetIdxW1Prod1 = -1; _genpartonJetIdxW1Prod2 = -1;
  _genpartonJetIdxB2 = -1; _genpartonJetIdxW2Prod1 = -1; _genpartonJetIdxW2Prod2 = -1;
  
  _recoTTBar  .clear();
  _recoTop1   .clear();
  _recoTop2   .clear();
  _recoW1     .clear();
  _recoW2     .clear();
  _recoB1     .clear();
  _recoW1Prod1.clear();
  _recoW1Prod2.clear();
  _recoB2     .clear();
  _recoW2Prod1.clear();
  _recoW2Prod2.clear();

  _recoJetIdxB1     .clear();
  _recoJetIdxW1Prod1.clear();
  _recoJetIdxW1Prod2.clear();
  _recoJetIdxB2     .clear();
  _recoJetIdxW2Prod1.clear();
  _recoJetIdxW2Prod2.clear();
  
  _fitTTBar    .clear();
  _fitTop1   .clear();
  _fitTop2   .clear();
  _fitW1     .clear();
  _fitW2     .clear();
  _fitB1     .clear();
  _fitW1Prod1.clear();
  _fitW1Prod2.clear();
  _fitB2     .clear();
  _fitW2Prod1.clear();
  _fitW2Prod2.clear();

  _combinationType.clear();

  _fitProb .clear();
  _fitChi2 .clear();
  _fitSigMT.clear();
}
