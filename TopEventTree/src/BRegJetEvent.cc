/*
 * BRegJetEvent.cc
 *
 *  Created on: May 14, 2013
 *      Author: kirschen
 */

#include "TopMass/TopEventTree/interface/BRegJetEvent.h"
#include "TopMass/TopEventTree/interface/BRegJetEvent_LinkDef.h"

BRegJetEvent::BRegJetEvent() : TObject()
{
  init();
}

void BRegJetEvent::init()
{
//  fChargedHadron.clear();
  EtWeightedSigmaPhi.clear();
  EtWeightedSigmaEta.clear();
  jesTotUnc.clear();
  jetPtRaw .clear();
  jetArea.clear();
  nChargedPFConstituents.clear();
  leadingChargedConstPt.clear();
  nSV.clear();
  SVChi2.clear();	  
  SV3DLength.clear();	  
  SV3DLengthError.clear();
  SVMass.clear();
  SVPt.clear();
  Rho.clear();
  Rho25.clear();

  nSoftMuons          .clear();
  nSoftElectrons      .clear();
  SoftMuonPt          .clear();
  SoftMuonPtRel       .clear();
  SoftMuonRatioRel       .clear();
  SoftMuonDeltaR      .clear();
//  SoftMuonJet_idx     .clear();
  SoftElectronPt      .clear();
  SoftElectronPtRel   .clear();
  SoftElectronRatioRel   .clear();
  SoftElectronDeltaR  .clear();
//  SoftElectronJet_idx .clear();

  BRegResult.clear();
}
