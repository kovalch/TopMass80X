#ifndef SherpaGenEvent_h
#define SherpaGenEvent_h

#include "TObject.h"
#include "TLorentzVector.h"
#include <vector>

class SherpaGenEvent : public TObject {
public:

  SherpaGenEvent();
  void init();

  ClassDef(SherpaGenEvent,1);

  // SherpaGenEvent data
  
  std::vector<double> flavour;
  std::vector<double> dr;
  
  std::vector<TLorentzVector> parton;
  std::vector<TLorentzVector> genJet;
};

#endif
