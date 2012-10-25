#ifndef MASSANALYZER_H
#define MASSANALYZER_H

#include <iostream>
#include <map>
#include <utility>

#include "TString.h"
#include "TTree.h"

#include <boost/program_options.hpp>

namespace po = boost::program_options;

class MassAnalyzer {
public:
  MassAnalyzer(const TString& identifier, TTree* tree);
  virtual ~MassAnalyzer() { delete fTree_; }
  
  
  virtual void Analyze(const TString& cuts, int i, int j) = 0;
  virtual void Analyze(const TString& cuts, int i, int j, po::variables_map vm) = 0;

  std::pair<double, double> GetValue(TString whichValue) const;
  std::map<TString, std::pair<double, double> > GetValues() const;

protected:
  void SetValue(TString whichValue, double val, double valError);

  TString fIdentifier_;
  TTree* fTree_;

  std::map<TString, std::pair<double, double> > values_;
};

#endif
