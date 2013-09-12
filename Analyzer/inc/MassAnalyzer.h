#ifndef MASSANALYZER_H
#define MASSANALYZER_H

#include <assert.h>
#include <iostream>
#include <map>
#include <utility>

#include "TString.h"
#include "TTree.h"

class MassAnalyzer {
public:
  MassAnalyzer(const TString& identifier, TTree* tree);
  virtual ~MassAnalyzer() { delete fTree_; }
  
  
  virtual void Analyze(const TString& cuts, int i, int j) = 0;

  std::pair<double, double> GetValue(const TString& whichValue) const;
  std::map<TString, std::pair<double, double> > GetValues() const;

protected:
  void SetValue(const TString& whichValue, double val, double valError);

  TString fIdentifier_;
  TTree* fTree_;

  std::map<TString, std::pair<double, double> > values_;
};

#endif
