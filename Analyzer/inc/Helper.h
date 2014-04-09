#ifndef HELPER_H
#define HELPER_H

#include <string>
#include <vector>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "Rtypes.h" // needed for Color_t and kBlack

#include "TH1F.h"
#include "TF1.h"

class TH1F;

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
  enum methodID {kGenMatch, kMVA, kIdeogram, kIdeogramNew, kIdeogramMin, kRooFit, kMaxMethods};
  static int methodID();
  static int getCMSEnergy();
  std::vector<double> readParameters(const char *whichParameter);
  std::vector<std::string> readParametersString(const char *whichParameter);
};

namespace HelperFunctions {
  // "best of" of Matthias' Usercode: https://github.com/mschrode/Tools
  TH1* createRatioPlot(const TH1 *h1, const TH1 *h2, const std::string &yTitle);
  std::string cleanedName(std::string toBeCleaned);
  void findYRange(const TH1 *h, double& min, double& max);
  void setCommonYRange(std::vector <TH1 *> histos, double RelTopOffset=0, double logMin =-1.);
  bool fitCoreWidth(const TH1* hist, double nSig, TF1* &gauss, double &width, double &widthErr, double &rms, double &rmsErr);
  bool equidistLogBins(std::vector<double>& binEdges, double min, double max, bool logarithm=true);
  std::string addProperArrayIndex(std::string inputexpression, std::string arrayIndex="");
}

#endif
