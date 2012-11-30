#include "TopMass/TopEventTree/interface/TopEvent.h"

void TopEvent::init() {
  run = -1; lumiBlock = -1; event = -1;
  
  genparton.clear();
  genpartonJetIdx.clear();
  
  jet.clear();
  jetCharge.clear();
  jetFlavour.clear();
  jetCSV.clear();
  
  reco.clear();
  recoJetIdx.clear();
  
  fit.clear();
  fitProb.clear();
  fitChi2.clear();
}

