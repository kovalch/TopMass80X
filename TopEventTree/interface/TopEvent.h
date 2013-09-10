#ifndef TopEvent_h
#define TopEvent_h

#include "TObject.h"
#include "TLorentzVector.h"
#include <vector>

class TopEvent : public TObject {
public:

  TopEvent();
  void init();

  ClassDef(TopEvent,6);

  // TopEvent data
  
  unsigned int run;
  unsigned int lumiBlock;
  unsigned int event;

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
