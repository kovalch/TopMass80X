/*
 * WeightEvent.cc
 *
 *  Created on: Feb 7, 2013
 *      Author: eschliec
 */

#include "TopMass/TopEventTree/interface/WeightEvent.h"

WeightEvent::WeightEvent()
{
  init();
}

void WeightEvent::init()
{
  mcWeight = 0.;

  puWeight     = -1.;
  puWeightUp   = -1.;
  puWeightDown = -1.;
  nVertex = -1;
  nPU.clear();
  nPUTrue = -1.;

  bTagEffJetWeight.clear();
  lTagEffJetWeight.clear();
  bTagEffWeight = -1.;
  lTagEffWeight = -1.;

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
