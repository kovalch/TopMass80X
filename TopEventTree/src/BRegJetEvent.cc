/*
 * BRegJetEvent.cc
 *
 *  Created on: May 14, 2013
 *      Author: kirschen
 */

#include "TopMass/TopEventTree/interface/BRegJetEvent.h"

BRegJetEvent::BRegJetEvent()
{
  init();
}

void BRegJetEvent::init()
{
  fChargedHadron.clear();
  EtWeightedSigmaPhi.clear();
  EtWeightedSigmaEta.clear();
  jesTotUnc.clear();
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

}
