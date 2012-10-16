#ifndef HELPER_H
#define HELPER_H

#include "TH2F.h"
#include "TStyle.h"
#include "TPaveLabel.h"

class Helper {
  private:
    int fBins;
    TString fBinning;
    std::vector<float> vBinning;
    void DrawLabel(TString text, const double x1, const double y1, const double x2, Color_t color);

  public:
    Helper(int bins, TString binning, std::vector<float> v)
      : fBins(bins), fBinning(binning), vBinning(v) {};
    Helper(int bins) : fBins(bins) {};
    ~Helper();
    
    TH1F* GetH1(TString title);
    TH2F* GetH2(TString title);
    void SetTDRStyle();
    void DrawCMSPrel();
    void DrawCMSPrelElectron();
    void DrawCMSPrelMuon();
    void DrawCMSSim();
};

#endif
