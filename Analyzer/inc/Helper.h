#ifndef HELPER_H
#define HELPER_H

#include "TH2F.h"
#include "TString.h"

class Helper {
  private:
    int fBins;
    void DrawLabel(TString text, const double x1, const double y1, const double x2, Color_t color = kBlack);

  public:
    Helper(int bins) : fBins(bins) {};
    ~Helper();
    
    TH2F* GetH2(TString title);
    void SetTDRStyle();
    void DrawCMSPrel();
    void DrawCMSSim();
};

#endif
