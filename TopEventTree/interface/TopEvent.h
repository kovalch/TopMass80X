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

  ClassDef(TopEvent,2);

  // TopEvent data
  
  int run;
  int lumiBlock;
  int event;
  
  int decayChannel;
  
  TLorentzVector genpartonTTBar;
  TLorentzVector genpartonHadTop;
  TLorentzVector genpartonLepTop;
  TLorentzVector genpartonHadW;
  TLorentzVector genpartonLepW;
  TLorentzVector genpartonHadB;
  TLorentzVector genpartonLightQ;
  TLorentzVector genpartonLightQBar;
  TLorentzVector genpartonLepB;
  TLorentzVector genpartonLepton;
  TLorentzVector genpartonNeutrino;
  
  TLorentzVector *genpartonTop1;
  TLorentzVector *genpartonTop2;
  TLorentzVector *genpartonW1;
  TLorentzVector *genpartonW2;
  TLorentzVector *genpartonB1;
  TLorentzVector *genpartonLightQ1;
  TLorentzVector *genpartonLightQBar1;
  TLorentzVector *genpartonB2;
  TLorentzVector *genpartonLightQ2;
  TLorentzVector *genpartonLightQBar2;
  
  int genpartonJetIdxHadB;
  int genpartonJetIdxLightQ;
  int genpartonJetIdxLightQBar;
  int genpartonJetIdxLepB;
  int genpartonJetIdxLepton;
  int genpartonJetIdxNeutrino;
  
  int genpartonJetIdxB1;
  int genpartonJetIdxLightQ1;
  int genpartonJetIdxLightQBar1;
  int genpartonJetIdxB2;
  int genpartonJetIdxLightQ2;
  int genpartonJetIdxLightQBar2;
  
  std::vector<TLorentzVector> recoTTBar;
  std::vector<TLorentzVector> recoHadTop;
  std::vector<TLorentzVector> recoLepTop;
  std::vector<TLorentzVector> recoHadW;
  std::vector<TLorentzVector> recoLepW;
  std::vector<TLorentzVector> recoHadB;
  std::vector<TLorentzVector> recoLightQ;
  std::vector<TLorentzVector> recoLightQBar;
  std::vector<TLorentzVector> recoLepB;
  std::vector<TLorentzVector> recoLepton;
  std::vector<TLorentzVector> recoNeutrino;

  std::vector<TLorentzVector> *recoTop1;
  std::vector<TLorentzVector> *recoTop2;
  std::vector<TLorentzVector> *recoW1;
  std::vector<TLorentzVector> *recoW2;
  std::vector<TLorentzVector> *recoB1;
  std::vector<TLorentzVector> *recoLightQ1;
  std::vector<TLorentzVector> *recoLightQBar1;
  std::vector<TLorentzVector> *recoB2;
  std::vector<TLorentzVector> *recoLightQ2;
  std::vector<TLorentzVector> *recoLightQBar2;

  std::vector<int> recoJetIdxHadB;
  std::vector<int> recoJetIdxLightQ;
  std::vector<int> recoJetIdxLightQBar;
  std::vector<int> recoJetIdxLepB;
  std::vector<int> recoJetIdxLepton;
  std::vector<int> recoJetIdxNeutrino;
  
  std::vector<int> *recoJetIdxB1;
  std::vector<int> *recoJetIdxLightQ1;
  std::vector<int> *recoJetIdxLightQBar1;
  std::vector<int> *recoJetIdxB2;
  std::vector<int> *recoJetIdxLightQ2;
  std::vector<int> *recoJetIdxLightQBar2;

  std::vector<TLorentzVector> fitTTBar;
  std::vector<TLorentzVector> fitHadTop;
  std::vector<TLorentzVector> fitLepTop;
  std::vector<TLorentzVector> fitHadW;
  std::vector<TLorentzVector> fitLepW;
  std::vector<TLorentzVector> fitHadB;
  std::vector<TLorentzVector> fitLightQ;
  std::vector<TLorentzVector> fitLightQBar;
  std::vector<TLorentzVector> fitLepB;
  std::vector<TLorentzVector> fitLepton;
  std::vector<TLorentzVector> fitNeutrino;

  std::vector<TLorentzVector> *fitTop1;
  std::vector<TLorentzVector> *fitTop2;
  std::vector<TLorentzVector> *fitW1;
  std::vector<TLorentzVector> *fitW2;
  std::vector<TLorentzVector> *fitB1;
  std::vector<TLorentzVector> *fitLightQ1;
  std::vector<TLorentzVector> *fitLightQBar1;
  std::vector<TLorentzVector> *fitB2;
  std::vector<TLorentzVector> *fitLightQ2;
  std::vector<TLorentzVector> *fitLightQBar2;

  std::vector<int> combinationType;

  std::vector<double> fitProb;
  std::vector<double> fitChi2;
  std::vector<double> fitSigMT;


private:
};

#endif
