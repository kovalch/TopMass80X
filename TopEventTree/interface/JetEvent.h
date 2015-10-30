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

  unsigned int nJet() const { return jet.size(); }

  ClassDef(JetEvent,3);

  // JetEvent data

  std::vector<TLorentzVector> jet;
  std::vector<TLorentzVector> genJet;
  std::vector<TLorentzVector> genParton;
  std::vector<TLorentzVector> alternativeJet;
  TLorentzVector met;

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
  std::vector<int>    nSV;
  std::vector<double> SVChi2;
  std::vector<double> SV3DLength;
  std::vector<double> SV3DLengthError;
  std::vector<TLorentzVector> SVMomentum;
  std::vector<double> gluonTag;

  std::vector<double> jerSF;
  std::vector<double> jesSF;
  std::vector<double> totalSF;
  std::vector<double> breg;

  double sumEt;

  std::vector<TVector2> pull;
  std::vector<TVector2> pullCharged;
  std::vector<TVector2> genPull;
  std::vector<TVector2> genPullCharged;


private:
};

#endif /* JETEVENT_H_ */
