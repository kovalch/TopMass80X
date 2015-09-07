#ifndef EventWeightBJES_h
#define EventWeightBJES_h

#include <memory>
#include <string>
#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TString.h"
#include "TVector2.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

class EventWeightBJES : public edm::EDProducer {

 public:
  explicit EventWeightBJES(const edm::ParameterSet&);
  ~EventWeightBJES();

 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob();

 private:
  edm::InputTag genJets_;
  edm::InputTag genParticles_;
  double nuDecayFractionSource_;
  double nuDecayFractionTarget_;

  edm::FileInPath fragSourceFile_;
  edm::FileInPath fragTargetFile_;

  TFile* sourceFile;
  TFile* targetFile;

  TH1F* sourceHist;
  TH1F* targetHist;

  /// histogram container for output
  std::map<std::string, TH1F*> hists;

};

#endif
