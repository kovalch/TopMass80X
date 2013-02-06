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
  jetCharge.clear();
  jetFlavour.clear();
  jetCSV.clear();
}
