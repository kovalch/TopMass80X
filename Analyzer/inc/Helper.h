#ifndef HELPER_H
#define HELPER_H

#include "TH2F.h"
#include "TStyle.h"
#include "TPaveLabel.h"

class Helper {
  private:
    TString fBinning;
    std::vector<float> vBinning;
    void DrawLabel(TString text, const double x1, const double y1, const double x2, Color_t color);

  public:
    Helper(TString binning, std::vector<float> v)
      : fBinning(binning), vBinning(v) {};
    Helper() {};
    ~Helper();
    
    TH1F* GetH1(TString title);
    void SetTDRStyle();
    void DrawCMSPrel();
    void DrawCMSPrelElectron();
    void DrawCMSPrelMuon();
    void DrawCMSSim();
};

#endif
