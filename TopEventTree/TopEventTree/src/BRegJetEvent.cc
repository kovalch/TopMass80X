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
  OneOfLeading2B.clear();
  RlbReco.clear();

  QGaxis1.clear();
  QGaxis2.clear();
  QGMult .clear();
  QGPtD  .clear();
  QGMLP  .clear();

  PUIddZ       .clear();
  PUIddRMean   .clear();
  PUIddr2Mean  .clear();
  PUIdfrac01   .clear();
  PUIdfrac02   .clear();
  PUIdfrac03   .clear();
  PUIdfrac04   .clear();
  PUIdfrac05   .clear();
  PUIdbeta     .clear();
  PUIdbetaStar .clear();
  PUIdptD      .clear();

  PUJetIdMVA   .clear();
  PUJetidflag  .clear();

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

  jetPtCorr     .clear();
  jetEta        .clear();
  jetMt         .clear();
  genJetPt      .clear();
  genPartonPt   .clear();


  BRegResult.clear();
  BRegGBRTrainResult.clear();
  BRegProducerResult.clear();
}
