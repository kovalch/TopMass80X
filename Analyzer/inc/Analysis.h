#include <map>
#include <string>
#include <vector>

#include "TTree.h"

class TH1F;
class MassAnalyzer;
class RandomSubsetCreator;

class Analysis {
 private:

  const std::string samplePath_, fIdentifier_, fMethod_;
  std::string fBinning_;
  const std::vector<float> vBinning_;
  const int fChannelID_, fMethodID_;
  MassAnalyzer* fAnalyzer_;
  RandomSubsetCreator* fCreator_;

  TTree* fTree_;

  std::map<std::string, TH1F*> histograms_;

  void SetH1(std::string histName, TH1F* hist);
  void CreateHisto(std::string name);

 public:
  Analysis(const std::vector<float>& v);
  ~Analysis();
    
  void Analyze();

  TH1F* GetH1(std::string histName);
  const std::map<std::string, TH1F*> GetH1s() const;

  std::string GetIdentifier();

};
