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
  leptonFlavour = 0;
  
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
  
  fitTTBar  .clear();
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

void TopEvent::shrink(unsigned int maxSize)
{
  if(recoTTBar.size() > maxSize){
    recoTTBar  .resize(maxSize);
    recoTop1   .resize(maxSize);
    recoTop2   .resize(maxSize);
    recoW1     .resize(maxSize);
    recoW2     .resize(maxSize);
    recoB1     .resize(maxSize);
    recoW1Prod1.resize(maxSize);
    recoW1Prod2.resize(maxSize);
    recoB2     .resize(maxSize);
    recoW2Prod1.resize(maxSize);
    recoW2Prod2.resize(maxSize);

    recoJetIdxB1     .resize(maxSize);
    recoJetIdxW1Prod1.resize(maxSize);
    recoJetIdxW1Prod2.resize(maxSize);
    recoJetIdxB2     .resize(maxSize);
    recoJetIdxW2Prod1.resize(maxSize);
    recoJetIdxW2Prod2.resize(maxSize);

    fitTTBar  .resize(maxSize);
    fitTop1   .resize(maxSize);
    fitTop2   .resize(maxSize);
    fitW1     .resize(maxSize);
    fitW2     .resize(maxSize);
    fitB1     .resize(maxSize);
    fitW1Prod1.resize(maxSize);
    fitW1Prod2.resize(maxSize);
    fitB2     .resize(maxSize);
    fitW2Prod1.resize(maxSize);
    fitW2Prod2.resize(maxSize);

    combinationType.resize(maxSize);

    fitProb.resize(maxSize);
    fitChi2.resize(maxSize);
    // treat sigMT specially as it is not used in all-jets channel up to now
    if(fitSigMT.size() > maxSize) fitSigMT.resize(maxSize);
  }
}
