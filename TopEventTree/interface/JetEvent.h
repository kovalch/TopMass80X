/*
 * JetEvent.h
 *
 *  Created on: Feb 6, 2013
 *      Author: eschliec
 */

#ifndef JETEVENT_H_
#define JETEVENT_H_

#include "TObject.h"
#include "TLorentzVector.h"
#include <vector>


class JetEvent : public TObject {
public:

  JetEvent();
  void init();

  unsigned int nJet() { return jet.size(); }

  ClassDef(JetEvent,1);

  // JetEvent data

  std::vector<TLorentzVector> jet;

  std::vector<int> nConstituents;
  std::vector<int> nChargedHadrons;
  std::vector<int> nNeutralHadrons;
  std::vector<int> nElectrons;
  std::vector<int> nPhotons;
  std::vector<int> nMuons;
  std::vector<double> fChargedHadron;
  std::vector<double> fNeutralHadron;
  std::vector<double> fElectron;
  std::vector<double> fPhoton;
  std::vector<double> fMuon;

  std::vector<double> charge;
  std::vector<double> flavour;
  std::vector<double> bTagCSV;
  std::vector<double> gluonTag;

  std::vector<double> jerSF;
  std::vector<double> jesSF;
  std::vector<double> totalSF;

  std::vector<TVector2> pull;
  std::vector<TVector2> pullCharged;


private:
};

#endif /* JETEVENT_H_ */
