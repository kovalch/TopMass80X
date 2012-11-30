#ifndef TopEvent_h
#define TopEvent_h

#include "TObject.h"
#include "TLorentzVector.h"
#include <vector>


class TopEvent : public TObject {
public:

  enum { TTBar, HadTop, LepTop, HadW, LepW, HadB, LightQ, LightQBar, LepB, Lepton, Neutrino };
  enum { TTBar0, Top1, Top2, W1, W2, B1, LightQ1, LightQBar1, B2, LightQ2, LightQBar2 };
  
  TopEvent() {};
  void init();

  double getFitTopMass(int permutation) { return fit[permutation][HadTop].M(); } // TODO crash after few entries
  
  ClassDef(TopEvent,1);

  // Event data
  
  int run;
  int lumiBlock;
  int event;
  
  std::vector<TLorentzVector> genparton;
  std::vector<int> genpartonJetIdx;
  
  std::vector<TLorentzVector> jet;
  std::vector<double> jetCharge;
  std::vector<double> jetFlavour;
  std::vector<double> jetCSV;
  
  std::vector<std::vector<TLorentzVector> > reco;
  std::vector<std::vector<int> > recoJetIdx;
  
  std::vector<std::vector<TLorentzVector> > fit;
  std::vector<double> fitProb;
  std::vector<double> fitChi2;


private:
};

#endif
