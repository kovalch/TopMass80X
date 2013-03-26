#ifndef HELPER_H
#define HELPER_H

#include "TH1F.h"
#include "TString.h"

#include <vector>

class Helper {
private:
  TString fBinning;
  std::vector<float> vBinning;
  void init();
  static int channelIDFromString(std::string channel);
  int energyFromString(std::string sample);
  void DrawLabel(TString text, const double x1, const double y1, const double x2, Color_t color = kBlack);
  int channelID_;
  int energy_;

public:
  Helper() {init();};
  Helper(TString binning, std::vector<float> v) : fBinning(binning), vBinning(v) {init();};
  ~Helper() {};

  TH1F* GetH1(TString title);
  void SetTDRStyle();
  void DrawCMS(int channelID = -1, int energy = -1);
  void DrawCMSSim(int energy = -1);

  enum channelID {kAllJets, kMuonJets, kElectronJets, kLeptonJets, kMaxChannels};
  static int channelID();
};

#endif
