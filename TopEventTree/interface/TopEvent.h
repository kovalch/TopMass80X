#ifndef TopEvent_h2Prod11Prod1
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

  ClassDef(TopEvent,3);

  // TopEvent data
  
  int& run      () {return _run      ;}
  int& lumiBlock() {return _lumiBlock;}
  int& event    () {return _event    ;}
  
  int& decayChannel() {return _decayChannel;}
  
  TLorentzVector& genpartonTTBar    () {return _genpartonTTBar  ;}
  TLorentzVector& genpartonHadTop   () {return _genpartonTop1   ;}
  TLorentzVector& genpartonLepTop   () {return _genpartonTop2   ;}
  TLorentzVector& genpartonHadW     () {return _genpartonW1     ;}
  TLorentzVector& genpartonLepW     () {return _genpartonW2     ;}
  TLorentzVector& genpartonHadB     () {return _genpartonB1     ;}
  TLorentzVector& genpartonLightQ   () {return _genpartonW1Prod1;}
  TLorentzVector& genpartonLightQBar() {return _genpartonW1Prod2;}
  TLorentzVector& genpartonLepB     () {return _genpartonB2     ;}
  TLorentzVector& genpartonLepton   () {return _genpartonW2Prod1;}
  TLorentzVector& genpartonNeutrino () {return _genpartonW2Prod2;}

  TLorentzVector& genpartonTop1      () {return _genpartonTop1   ;}
  TLorentzVector& genpartonTop2      () {return _genpartonTop2   ;}
  TLorentzVector& genpartonW1        () {return _genpartonW1     ;}
  TLorentzVector& genpartonW2        () {return _genpartonW2     ;}
  TLorentzVector& genpartonB1        () {return _genpartonB1     ;}
  TLorentzVector& genpartonLightQ1   () {return _genpartonW1Prod1;}
  TLorentzVector& genpartonLightQBar1() {return _genpartonW1Prod2;}
  TLorentzVector& genpartonB2        () {return _genpartonB2     ;}
  TLorentzVector& genpartonLightQ2   () {return _genpartonW2Prod1;}
  TLorentzVector& genpartonLightQBar2() {return _genpartonW2Prod2;}
  
  TLorentzVector& genpartonLepton1  () {return _genpartonW1Prod1;}
  TLorentzVector& genpartonNeutrino1() {return _genpartonW1Prod2;}
  TLorentzVector& genpartonLepton2  () {return _genpartonW2Prod1;}
  TLorentzVector& genpartonNeutrino2() {return _genpartonW2Prod2;}

  int& genpartonJetIdxHadB     () {return _genpartonJetIdxB1     ;}
  int& genpartonJetIdxLightQ   () {return _genpartonJetIdxW1Prod1;}
  int& genpartonJetIdxLightQBar() {return _genpartonJetIdxW1Prod2;}
  int& genpartonJetIdxLepB     () {return _genpartonJetIdxB2     ;}
  int& genpartonJetIdxLepton   () {return _genpartonJetIdxW2Prod1;}
  int& genpartonJetIdxNeutrino () {return _genpartonJetIdxW2Prod2;}
  
  int& genpartonJetIdxB1        () {return _genpartonJetIdxB1     ;}
  int& genpartonJetIdxLightQ1   () {return _genpartonJetIdxW1Prod1;}
  int& genpartonJetIdxLightQBar1() {return _genpartonJetIdxW1Prod2;}
  int& genpartonJetIdxB2        () {return _genpartonJetIdxB2     ;}
  int& genpartonJetIdxLightQ2   () {return _genpartonJetIdxW2Prod1;}
  int& genpartonJetIdxLightQBar2() {return _genpartonJetIdxW2Prod2;}
  
  std::vector<TLorentzVector>& recoTTBar    () {return _recoTTBar  ;}
  std::vector<TLorentzVector>& recoHadTop   () {return _recoTop1   ;}
  std::vector<TLorentzVector>& recoLepTop   () {return _recoTop2   ;}
  std::vector<TLorentzVector>& recoHadW     () {return _recoW1     ;}
  std::vector<TLorentzVector>& recoLepW     () {return _recoW2     ;}
  std::vector<TLorentzVector>& recoHadB     () {return _recoB1     ;}
  std::vector<TLorentzVector>& recoLightQ   () {return _recoW1Prod1;}
  std::vector<TLorentzVector>& recoLightQBar() {return _recoW1Prod2;}
  std::vector<TLorentzVector>& recoLepB     () {return _recoB2     ;}
  std::vector<TLorentzVector>& recoLepton   () {return _recoW2Prod1;}
  std::vector<TLorentzVector>& recoNeutrino () {return _recoW2Prod2;}

  std::vector<TLorentzVector>& recoTop1      () {return _recoTop1   ;}
  std::vector<TLorentzVector>& recoTop2      () {return _recoTop2   ;}
  std::vector<TLorentzVector>& recoW1        () {return _recoW1     ;}
  std::vector<TLorentzVector>& recoW2        () {return _recoW2     ;}
  std::vector<TLorentzVector>& recoB1        () {return _recoB1     ;}
  std::vector<TLorentzVector>& recoLightQ1   () {return _recoW1Prod1;}
  std::vector<TLorentzVector>& recoLightQBar1() {return _recoW1Prod2;}
  std::vector<TLorentzVector>& recoB2        () {return _recoB2     ;}
  std::vector<TLorentzVector>& recoLightQ2   () {return _recoW2Prod1;}
  std::vector<TLorentzVector>& recoLightQBar2() {return _recoW2Prod2;}

  std::vector<int>& recoJetIdxHadB     () {return _recoJetIdxB1     ;}
  std::vector<int>& recoJetIdxLightQ   () {return _recoJetIdxW1Prod1;}
  std::vector<int>& recoJetIdxLightQBar() {return _recoJetIdxW1Prod2;}
  std::vector<int>& recoJetIdxLepB     () {return _recoJetIdxB2     ;}
  std::vector<int>& recoJetIdxLepton   () {return _recoJetIdxW2Prod1;}
  std::vector<int>& recoJetIdxNeutrino () {return _recoJetIdxW2Prod2;}
  
  std::vector<int>& recoJetIdxB1        () {return _recoJetIdxB1     ;}
  std::vector<int>& recoJetIdxLightQ1   () {return _recoJetIdxW1Prod1;}
  std::vector<int>& recoJetIdxLightQBar1() {return _recoJetIdxW1Prod2;}
  std::vector<int>& recoJetIdxB2        () {return _recoJetIdxB2     ;}
  std::vector<int>& recoJetIdxLightQ2   () {return _recoJetIdxW2Prod1;}
  std::vector<int>& recoJetIdxLightQBar2() {return _recoJetIdxW2Prod2;}

  std::vector<TLorentzVector>& fitTTBar    () {return _fitTTBar  ;}
  std::vector<TLorentzVector>& fitHadTop   () {return _fitTop1   ;}
  std::vector<TLorentzVector>& fitLepTop   () {return _fitTop2   ;}
  std::vector<TLorentzVector>& fitHadW     () {return _fitW1     ;}
  std::vector<TLorentzVector>& fitLepW     () {return _fitW2     ;}
  std::vector<TLorentzVector>& fitHadB     () {return _fitB1     ;}
  std::vector<TLorentzVector>& fitLightQ   () {return _fitW1Prod1;}
  std::vector<TLorentzVector>& fitLightQBar() {return _fitW1Prod2;}
  std::vector<TLorentzVector>& fitLepB     () {return _fitB2     ;}
  std::vector<TLorentzVector>& fitLepton   () {return _fitW2Prod1;}
  std::vector<TLorentzVector>& fitNeutrino () {return _fitW2Prod2;}

  std::vector<TLorentzVector>& fitTop1      () {return _fitTop1   ;}
  std::vector<TLorentzVector>& fitTop2      () {return _fitTop2   ;}
  std::vector<TLorentzVector>& fitW1        () {return _fitW1     ;}
  std::vector<TLorentzVector>& fitW2        () {return _fitW2     ;}
  std::vector<TLorentzVector>& fitB1        () {return _fitB1     ;}
  std::vector<TLorentzVector>& fitLightQ1   () {return _fitW1Prod1;}
  std::vector<TLorentzVector>& fitLightQBar1() {return _fitW1Prod2;}
  std::vector<TLorentzVector>& fitB2        () {return _fitB2     ;}
  std::vector<TLorentzVector>& fitLightQ2   () {return _fitW2Prod1;}
  std::vector<TLorentzVector>& fitLightQBar2() {return _fitW2Prod2;}

  std::vector<int>& combinationType() {return _combinationType;}

  std::vector<double>& fitProb () {return _fitProb;}
  std::vector<double>& fitChi2 () {return _fitChi2;}
  std::vector<double>& fitSigMT() {return _fitSigMT;}


private:

  int _run;
  int _lumiBlock;
  int _event;

  int _decayChannel;

  TLorentzVector _genpartonTTBar;
  TLorentzVector _genpartonTop1;
  TLorentzVector _genpartonTop2;
  TLorentzVector _genpartonW1;
  TLorentzVector _genpartonW2;
  TLorentzVector _genpartonB1;
  TLorentzVector _genpartonW1Prod1;
  TLorentzVector _genpartonW1Prod2;
  TLorentzVector _genpartonB2;
  TLorentzVector _genpartonW2Prod1;
  TLorentzVector _genpartonW2Prod2;

  int _genpartonJetIdxB1;
  int _genpartonJetIdxW1Prod1;
  int _genpartonJetIdxW1Prod2;
  int _genpartonJetIdxB2;
  int _genpartonJetIdxW2Prod1;
  int _genpartonJetIdxW2Prod2;

  std::vector<TLorentzVector> _recoTTBar;
  std::vector<TLorentzVector> _recoTop1;
  std::vector<TLorentzVector> _recoTop2;
  std::vector<TLorentzVector> _recoW1;
  std::vector<TLorentzVector> _recoW2;
  std::vector<TLorentzVector> _recoB1;
  std::vector<TLorentzVector> _recoW1Prod1;
  std::vector<TLorentzVector> _recoW1Prod2;
  std::vector<TLorentzVector> _recoB2;
  std::vector<TLorentzVector> _recoW2Prod1;
  std::vector<TLorentzVector> _recoW2Prod2;

  std::vector<int> _recoJetIdxB1;
  std::vector<int> _recoJetIdxW1Prod1;
  std::vector<int> _recoJetIdxW1Prod2;
  std::vector<int> _recoJetIdxB2;
  std::vector<int> _recoJetIdxW2Prod1;
  std::vector<int> _recoJetIdxW2Prod2;

  std::vector<TLorentzVector> _fitTTBar;
  std::vector<TLorentzVector> _fitTop1;
  std::vector<TLorentzVector> _fitTop2;
  std::vector<TLorentzVector> _fitW1;
  std::vector<TLorentzVector> _fitW2;
  std::vector<TLorentzVector> _fitB1;
  std::vector<TLorentzVector> _fitW1Prod1;
  std::vector<TLorentzVector> _fitW1Prod2;
  std::vector<TLorentzVector> _fitB2;
  std::vector<TLorentzVector> _fitW2Prod1;
  std::vector<TLorentzVector> _fitW2Prod2;

  std::vector<int> _combinationType;

  std::vector<double> _fitProb;
  std::vector<double> _fitChi2;
  std::vector<double> _fitSigMT;

};

#endif
