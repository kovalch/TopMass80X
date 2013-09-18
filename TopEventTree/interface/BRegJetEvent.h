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

  unsigned int nJet() { return EtWeightedSigmaPhi.size(); }

  ClassDef(BRegJetEvent,2);

  // BRegJetEvent data

//  std::vector<double> fChargedHadron;
  std::vector<double> EtWeightedSigmaPhi;
  std::vector<double> EtWeightedSigmaEta;
  std::vector<double> jesTotUnc;
  std::vector<double> jetPtRaw;
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

  std::vector<int>    nSoftMuons;
  std::vector<int>    nSoftElectrons;
  std::vector<double> SoftMuonPt;
  std::vector<double> SoftMuonPtRel;         // transverse momentum wrt. the jet axis
  std::vector<double> SoftMuonRatioRel;         // momentum parallel to jet axis over jet energy
  std::vector<double> SoftMuonDeltaR;// (pseudo)angular distance to jet axis
//  std::vector<double> SoftMuonJet_idx;
  std::vector<double> SoftElectronPt;
  std::vector<double> SoftElectronPtRel;             // transverse momentum wrt. the jet axis
  std::vector<double> SoftElectronRatioRel;             // momentum parallel to jet axis over jet energy
  std::vector<double> SoftElectronDeltaR;// (pseudo)angular distance to jet axis
//  std::vector<double> SoftElectronJet_idx;


  std::vector<double> jetPtCorr;
  std::vector<double> jetEta;
  std::vector<double> jetMt;
  std::vector<double> genJetPt;
  std::vector<double> genPartonPt;


  std::vector<double> BRegResult;
  std::vector<double> BRegGBRTrainResult;


private:
};

#endif /* BREGJETEVENT_H_ */
