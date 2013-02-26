/*
 * JetEvent.cc
 *
 *  Created on: Feb 6, 2013
 *      Author: eschliec
 */

#include "TopMass/TopEventTree/interface/JetEvent.h"

JetEvent::JetEvent()
{
  init();
}

void JetEvent::init()
{
  jet.clear();
  genJet.clear();
  genParton.clear();

  nConstituents.clear();
  nChargedHadrons.clear();
  nNeutralHadrons.clear();
  nElectrons.clear();
  nPhotons.clear();
  nMuons.clear();
  fChargedHadron.clear();
  fNeutralHadron.clear();
  fElectron.clear();
  fPhoton.clear();
  fMuon.clear();

  charge.clear();
  flavour.clear();
  bTagCSV.clear();
  gluonTag.clear();

  jerSF.clear();
  jesSF.clear();
  totalSF.clear();

  // CINT does not like vectors of TVector2s
  // But if you don't tell ROOT, what they are, it works :D
#ifndef __CINT__
  pull.clear();
  pullCharged.clear();
#endif
}
