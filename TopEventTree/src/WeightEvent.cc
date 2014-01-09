/*
 * WeightEvent.cc
 *
 *  Created on: Feb 7, 2013
 *      Author: eschliec
 */

#include "TopMass/TopEventTree/interface/WeightEvent.h"
#include "TopMass/TopEventTree/interface/WeightEvent_LinkDef.h"

WeightEvent::WeightEvent() : TObject()
{
  init();
}

void WeightEvent::init()
{
  mcWeight = 0.;
  brWeight = 1.;

  puWeight     = -1.;
  puWeightUp   = -1.;
  puWeightDown = -1.;
  nVertex = -1;
  nPU.clear();
  nPUTrue = -1.;

  bTagEffJetWeight.clear();
  lTagEffJetWeight.clear();
  bTagWeight              = -1.;
  bTagWeight_bTagSFUp     = -1.;
  bTagWeight_bTagSFDown   = -1.;
  bTagWeight_misTagSFUp   = -1.;
  bTagWeight_misTagSFDown = -1.;
  
  bJESWeight_fNuUp     = -1.;
  bJESWeight_fNuDown   = -1.;
  bJESWeight_frag      = -1.;
  bJESWeight_fragHard  = -1.;
  bJESWeight_fragSoft  = -1.;

  triggerWeight = -1.;

  muWeight = -1.;
  elWeight = -1.;

  pdfWeight.clear();
  x1 = -1.;
  x2 = -1.;
  Q  = -1.;
  id1 = -100;
  id2 = -100;

  combinedWeight = 1.;
}
