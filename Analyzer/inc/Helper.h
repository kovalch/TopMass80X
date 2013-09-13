#ifndef HELPER_H
#define HELPER_H

#include "TH1F.h"

#include <string>
#include <vector>

class Helper {
private:
  const std::string fBinning;
  const std::vector<float> vBinning;
  void init();
  static int channelIDFromString(const std::string& channel);
  static int methodIDFromString(const std::string& method);
  void DrawLabel(const std::string& text, const double x1, const double y1, const double x2, Color_t color = kBlack);
  int channelID_;
  int energy_;

public:
  Helper() {init();};
  Helper(const std::string& binning, const std::vector<float>& v) : fBinning(binning), vBinning(v) {init();};
  ~Helper() {};

  TH1F* GetH1(const std::string& title);
  void SetTDRStyle();
  void DrawCMS(int channelID = -1, int energy = -1);
  void DrawCMSSim(int energy = -1);

  enum channelID {kAllJets, kMuonJets, kElectronJets, kLeptonJets, kMaxChannels};
  static int channelID();
  enum methodID {kGenMatch, kMVA, kIdeogram, kIdeogramNew, kRooFit, kMaxMethods};
  static int methodID();
};

#endif
