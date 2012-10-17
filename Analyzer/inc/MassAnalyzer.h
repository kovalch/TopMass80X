#ifndef MASSANALYZER_H
#define MASSANALYZER_H

#include <iostream>

#include "TString.h"
#include "TTree.h"

#include <boost/program_options.hpp>

namespace po = boost::program_options;

class MassAnalyzer {
public:
  MassAnalyzer(const TString& identifier, TTree* tree);
  virtual ~MassAnalyzer() { delete fTree; }
  
  
  virtual void Analyze(const TString& cuts, int i, int j) = 0;
  virtual void Analyze(const TString& cuts, int i, int j, po::variables_map vm) = 0;
  
  double GetMass() const;
  double GetMassError() const;
  //double GetMassSigma() const;
  double GetJES() const;
  double GetJESError() const;

  double GetMassConstJES() const;
  double GetMassConstJESError() const;
  //double GetMassConstJESSigma() const;

  double GetFSig() const;
  double GetFSigError() const;
  double GetMassfSig() const;
  double GetMassfSigError() const;
  //double GetMassfSigSigma() const;
  double GetJESfSig() const;
  double GetJESfSigError() const;
  
protected:
  TString fIdentifier;
  TTree* fTree;
  double fEntries;
  double fMass;
  double fMassError;
  //double fMassSigma;
  double fJES;
  double fJESError;

  double fMassConstJES;
  double fMassConstJESError;
  //double fMassConstJESSigma;

  double ffSig;
  double ffSigError;
  double fMassfSig;
  double fMassfSigError;
  //double fMassfSigSigma;
  double fJESfSig;
  double fJESfSigError;
};

#endif
