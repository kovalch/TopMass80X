#ifndef MASSANALYZER_H
#define MASSANALYZER_H

#include <assert.h>
#include <iostream>
#include <map>
#include <string>
#include <utility>

#include "TTree.h"

class MassAnalyzer {
public:
  MassAnalyzer(const std::string& identifier, TTree* tree);
  virtual ~MassAnalyzer() { delete fTree_; }
  
  
  virtual void Analyze(const std::string& cuts, int i, int j) = 0;

  std::pair<double, double> GetValue(const std::string& whichValue) const;
  std::map<std::string, std::pair<double, double> > GetValues() const;

protected:
  void SetValue(const std::string& whichValue, double val, double valError);

  std::string fIdentifier_;
  TTree* fTree_;

  std::map<std::string, std::pair<double, double> > values_;
};

#endif
