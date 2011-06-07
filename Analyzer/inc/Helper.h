#ifndef HELPER_H
#define HELPER_H

#include "TH2F.h"
#include "TStyle.h"

class Helper {
  private:
    int fBins;

  public:
    Helper(int bins) : fBins(bins) {};
    ~Helper();
    
    TH2F* GetH2(TString title);
    void SetTDRStyle();
};

#endif
