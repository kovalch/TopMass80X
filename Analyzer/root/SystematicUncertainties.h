#ifndef SYSTEMATICUNCERTAINTIES_H
#define SYSTEMATICUNCERTAINTIES_H

#include <map>
#include <string>
#include <vector>

#include "TChain.h"

class SystematicUncertainties {
public:
  SystematicUncertainties();
  ~SystematicUncertainties() {}

private:

  struct ensemble {
    std::string file; // file selection for chain
    double size; // EFFECTIVE sample size
    std::map<std::string, std::pair<double, double>> values;
    TChain* chain;
    ensemble(std::string f = "temp", double s = 0, std::vector<std::pair<std::string, std::pair<double,double>>> v = std::vector<std::pair<std::string, std::pair<double,double>>>(0))
    : file(f), size(s), chain(0)
    {
      for(auto& var : v){
        values[var.first] = var.second;
      }
    }
  };

  struct comparison {
    std::string nominal;
    std::string up;
    std::string down;
    std::map<std::string, double> shifts;
    std::map<std::string, double> shiftUncs;
    bool correlated;
    bool active;
    comparison(std::string n = "", std::string u = "", std::string d = "", bool c = true, bool a = true)
    : nominal(n), up(u), down(d), correlated(c), active(a) {}
  };

  struct mergedcomparison {
    std::vector<std::string> comparisons;
    mergedcomparison()
    : comparisons(0) {}
    mergedcomparison(std::vector<std::string> c)
    : comparisons(c) {}
  };

  struct dataSample {
    std::map<std::string, ensemble> ensembles;
    std::map<std::string, comparison> comparisons;
    std::map<std::string, mergedcomparison> mergedcomparisons;
    std::vector<std::string> variables;
    std::string path;
    double peLumi;
    double crossSection;
  };

  void deriveSystematics();
  void fillLeptonJets();
  void fillAllJets();

  std::string getSelectionFromVariable(std::string& var);

  dataSample sample;

};

#endif /* SYSTEMATICUNCERTAINTIES_H */

