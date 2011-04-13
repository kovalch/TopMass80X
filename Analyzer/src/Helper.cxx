#include "Helper.h"

TH2F* Helper::GetH2(TString title) {
  TH2F* hHelper = new TH2F();
  hHelper->SetBins(fBins, 0, 3, fBins, 0, 3);
  hHelper->SetStats(false);
  hHelper->SetTitle(title);
  hHelper->SetXTitle("#theta^{decay}_{t}");
  hHelper->SetYTitle("#theta^{decay}_{W}");

  return hHelper;
}
