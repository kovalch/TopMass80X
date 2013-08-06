#ifndef TopEvent_h
#define TopEvent_h

#include "TObject.h"
#include "TLorentzVector.h"
#include <vector>

class TopEvent : public TObject {
public:

  enum semiLepParticles{ TTBar, HadTop, LepTop, HadW, LepW, HadB, LightQ, LightQBar, LepB, Lepton, Neutrino };
  enum fullHadParticles{ TTBar0, Top1, Top2, W1, W2, B1, LightQ1, LightQBar1, B2, LightQ2, LightQBar2 };
  
  TopEvent();
  void init();

  ClassDef(TopEvent,5);

  // TopEvent data
  
  int run;
  int lumiBlock;
  int event;

  int decayChannel;
  int leptonFlavour;

  TLorentzVector genpartonTTBar;
  TLorentzVector genpartonTop1;
  TLorentzVector genpartonTop2;
  TLorentzVector genpartonW1;
  TLorentzVector genpartonW2;
  TLorentzVector genpartonB1;
  TLorentzVector genpartonW1Prod1;
  TLorentzVector genpartonW1Prod2;
  TLorentzVector genpartonB2;
  TLorentzVector genpartonW2Prod1;
  TLorentzVector genpartonW2Prod2;

  int genpartonJetIdxB1;
  int genpartonJetIdxW1Prod1;
  int genpartonJetIdxW1Prod2;
  int genpartonJetIdxB2;
  int genpartonJetIdxW2Prod1;
  int genpartonJetIdxW2Prod2;

  std::vector<TLorentzVector> recoTTBar;
  std::vector<TLorentzVector> recoTop1;
  std::vector<TLorentzVector> recoTop2;
  std::vector<TLorentzVector> recoW1;
  std::vector<TLorentzVector> recoW2;
  std::vector<TLorentzVector> recoB1;
  std::vector<TLorentzVector> recoW1Prod1;
  std::vector<TLorentzVector> recoW1Prod2;
  std::vector<TLorentzVector> recoB2;
  std::vector<TLorentzVector> recoW2Prod1;
  std::vector<TLorentzVector> recoW2Prod2;

  std::vector<int> recoJetIdxB1;
  std::vector<int> recoJetIdxW1Prod1;
  std::vector<int> recoJetIdxW1Prod2;
  std::vector<int> recoJetIdxB2;
  std::vector<int> recoJetIdxW2Prod1;
  std::vector<int> recoJetIdxW2Prod2;

  std::vector<TLorentzVector> fitTTBar;
  std::vector<TLorentzVector> fitTop1;
  std::vector<TLorentzVector> fitTop2;
  std::vector<TLorentzVector> fitW1;
  std::vector<TLorentzVector> fitW2;
  std::vector<TLorentzVector> fitB1;
  std::vector<TLorentzVector> fitW1Prod1;
  std::vector<TLorentzVector> fitW1Prod2;
  std::vector<TLorentzVector> fitB2;
  std::vector<TLorentzVector> fitW2Prod1;
  std::vector<TLorentzVector> fitW2Prod2;

  std::vector<int> combinationType;

  std::vector<double> fitProb;
  std::vector<double> fitChi2;
  std::vector<double> fitSigMT;

};

#endif
