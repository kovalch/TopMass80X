#ifndef SemiLepLeptonAnalyzer_h
#define SemiLepLeptonAnalyzer_h

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"


/**
   \class   SemiLepLeptonAnalyzer.cc SemiLepLeptonAnalyzer.h "TopAnalysis/TopAnalyzer/plugins/SemiLepLeptonAnalyzer.h"

   \brief   Derived class to analyze the kinematics of leptons on reconstruction level (kinfit output) and generator level (leptons from genParticle collection)
   
   The structure keeps histograms for the kinematics of leptons. These histograms are filled either from the TtSemiLeptonicEvent class (reconstruction level) or from the specified genParticle collection. In the latter case the leading object is used. Furthermore, also the correlation plot between reconstructed and generated quantity is filled. 
*/

class SemiLepLeptonAnalyzer : public edm::EDAnalyzer {

 public:

  explicit SemiLepLeptonAnalyzer(const edm::ParameterSet&);
  ~SemiLepLeptonAnalyzer(){};
  
 private:

  // ---
  //    functions
  // ---  
  // default (framework)
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  // helper functions to get lepton from kinfit hypothesis
  const reco::Candidate* getLepton(const edm::Handle<TtSemiLeptonicEvent> semiLepEvt, const std::string hypoKey_);
  
  // ---
  //    inputs
  // ---
  // kinfit output
  edm::InputTag semiLepEvt_;
  // class key of kinfit hypothesis
  std::string hypoKey_;
  // reco Lepton collection 
  edm::InputTag recLeptons_;
  // genLepton collection 
  edm::InputTag genLeptons_;
  // output manager
  int verbose;
  // event weight
  edm::InputTag weight_;
  // create gen plots?
  bool genPlots_;
  // within visible cross section range?
  edm::InputTag ingenPS_;
  // create rec plots?
  bool recPlots_;
  // use lepton kinematics as before the kinematic fit?
  bool preLep;
  // event is within gen visible PS?
  bool inVisPS;

  // ---
  //    define Tree for event content
  // ---
  TTree * tree;
  bool useTree_;
  float weight;
  unsigned int runNumber, luminosityBlockNumber, eventNumber;
  float Q;
  int id1, id2;
  double x1, x2;
  float valueLepPtRec;
  float valueLepPtGen;
  float valueLepEtaRec;
  float valueLepEtaGen;
  float valueLepYRec;
  float valueLepYGen;

  // ---
  //    histos
  // ---
  TH1F* lepPtRec;
  TH1F* lepPtGen;

  TH1F* lepEtaRec;
  TH1F* lepEtaGen;

  TH1F* lepYRec;
  TH1F* lepYGen;

  TH2F* lepPt_;
  TH2F* lepEta_;
  TH2F* lepY_;

};

#endif
