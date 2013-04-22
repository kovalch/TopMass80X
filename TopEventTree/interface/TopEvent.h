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

  ClassDef(TopEvent,1);

  // TopEvent data
  
  int run;
  int lumiBlock;
  int event;
  
  int decayChannel;
  std::vector<int> genpartonJetIdx;
  std::vector<TLorentzVector> genparton;
  
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


private:
};

#endif
