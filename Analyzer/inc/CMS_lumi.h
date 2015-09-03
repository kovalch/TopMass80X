#ifndef CMSLUMI_H
#define CMSLUMI_H

#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TASImage.h"

//
// Global variables
//

class CMS_lumi {

public:
  CMS_lumi() {};
  void Draw_CMS_lumi( TPad* pad, int iPeriod=3, int iPosX=10 );

};

#endif
