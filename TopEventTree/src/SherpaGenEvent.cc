#include "TopMass/TopEventTree/interface/SherpaGenEvent.h"
#include "TopMass/TopEventTree/interface/SherpaGenEvent_LinkDef.h"

SherpaGenEvent::SherpaGenEvent() : TObject()
{
  init();
}

void SherpaGenEvent::init()
{
  flavour.clear();
  dr.clear();
  
  parton.clear();
  genJet.clear();
}

