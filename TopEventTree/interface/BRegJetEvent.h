/*
 * BRegJetEvent.h
 *
 *  Created on: May 14, 2013
 *      Author: kirschen
 */

#ifndef BREGJETEVENT_H_
#define BREGJETEVENT_H_

#include "TObject.h"
#include "TLorentzVector.h"
#include <vector>

class BRegJetEvent : public TObject {
public:

  BRegJetEvent();
  void init();

  unsigned int nJet() { return fChargedHadron.size(); }

  ClassDef(BRegJetEvent,1);

  // BRegJetEvent data

  std::vector<double> fChargedHadron;
  std::vector<double> EtWeightedSigmaPhi;
  std::vector<double> EtWeightedSigmaEta;
  std::vector<double> jesTotUnc;
  std::vector<double> jetArea;
  std::vector<int>    nChargedPFConstituents;
  std::vector<double> leadingChargedConstPt;
  std::vector<int>    nSV;
  std::vector<double> SVChi2;
  std::vector<double> SV3DLength;
  std::vector<double> SV3DLengthError;
  std::vector<double> SVMass;
  std::vector<double> SVPt;
  std::vector<double> Rho;
  std::vector<double> Rho25;


private:
};

#endif /* BREGJETEVENT_H_ */
