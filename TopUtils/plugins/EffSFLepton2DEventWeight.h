#ifndef EffSFMuonEventWeight_h
#define EffSFMuonEventWeight_h

#include <memory>
#include <string>
#include <iostream>

#include "TH1.h"
#include "TH2.h" 
#include "TFile.h"
#include "TString.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

/// This module returns an efficiency scale factor (SF)
/// which is put into the CMSSW event
/// as a double, which can be used as an event weight in the analyzers of interest.
/// cfg parameters:
/// particles_             particle collection
/// sysVar_                name of the systematic variation 
///                        "noSys"
///                        "combinedEffSFNormUpStat/Down": normalisation uncertainty using the statistical errors
///                        "combinedEffSFShapeUpEta/Down" : uncertainty to distort the shape (softens or tightens the dependence)
///                        "flatTriggerSF": uses meanTriggerEffSF_ as flat SF
///                        "combinedEffSFNormUpSys/Down": uses additionalFactorErr_ as flat global uncertainty
/// verbose_               set to 0 if no output on terminal is desired, 1 for moderate and 2 for detailed output
/// filename_              if not set to "", efficiencies are loaded from histos in filename_
/// additionalSystErr_      multiplies with factor ( + 0.5% for single muon triggers to be applyied on top of statistical errors)


class EffSFLepton2DEventWeight : public edm::EDProducer {

 public:
  explicit EffSFLepton2DEventWeight(const edm::ParameterSet&);
  ~EffSFLepton2DEventWeight();
  
 private:
//   virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob();

 private:
  edm::InputTag particles_;
  std::string sysVar_;
  int verbose_;
  edm::FileInPath filename_;
  std::string histoname_;
  double additionalSystErr_;
  double etaMaxVal_;
  
  /// helper values from 2D SFhisto
  double ptmin ;
  double ptmax ;
  double etamin;
  double etamax;
  
  /// histogram
  TH2F * hists2D;
  
  /// histogram container
  std::map<std::string, TH1F*> hists_;
  
  /// file with histos
  TFile * file_;
};

#endif
