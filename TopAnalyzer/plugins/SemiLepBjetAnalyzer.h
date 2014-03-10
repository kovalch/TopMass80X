#ifndef SemiLepBjetAnalyzer_h
#define SemiLepBjetAnalyzer_h

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"


/**
   \class   SemiLepBjetAnalyzer.cc SemiLepBjetAnalyzer.h "TopAnalysis/TopAnalyzer/interface/SemiLepBjetAnalyzer.h"

   \brief   Derived class to analyze the kinematics of b-jets on reconstruction level (kinfit output) and generator level (genjets from)
   
   The structure keeps histograms for the kinematics of bjets. These histograms are filled either from the TtSemiLeptonicEvent class (reconstruction level) or from the specified genjet collection. In the latter case the identification as done in the GenLevelBJetProducer (from TopAnalysis/TopUtils) is used. Furthermore, also the correlation plot between reconstructed and generated quantity is filled. The lepton charge is used to identify the b and bbar jets from the kinematic fit output. Both (b and bbar) are always drawn in the same histogram.   
*/

class SemiLepBjetAnalyzer : public edm::EDAnalyzer {

 public:

  explicit SemiLepBjetAnalyzer(const edm::ParameterSet&);
  ~SemiLepBjetAnalyzer(){};
  
 private:

  // ---
  //    functions
  // ---  
  // default (framework)
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  // helper functions to get b and bbar jet 
  // a) from kinFit hypothesis
  std::pair<const reco::Candidate*,const reco::Candidate*> getbJets(const edm::Handle<TtSemiLeptonicEvent> semiLepEvt, const std::string hypoKey_);
  // b) from kinFit assignment without kinematic shift
  std::pair<const reco::Candidate*,const reco::Candidate*> getPreBJets(const edm::Handle<TtSemiLeptonicEvent> semiLepEvt, const std::string hypoKey_, const edm::Handle<std::vector<pat::Jet> > jets);
  // c) from genJet collection
  const reco::GenJet* getJetFromCollection(const reco::GenJetCollection& genJets, int IX);
  // helper function to get the jet-parton matching result for the KinFit reco jet assignment
  double checkPartonAssignment(const edm::Handle<TtSemiLeptonicEvent> semiLepEvt, int maxNJets);

  // ---
  //    inputs
  // ---
  // kinfit output
  edm::InputTag semiLepEvt_;
  // class key of kinfit hypothesis
  std::string hypoKey_;
  // genJet collection 
  edm::InputTag genJets_;
  // genLepton collection 
  edm::InputTag genLeptons_;
  // pure b gen jet collection?
  bool bJetCollection_;
  // recoJet collection 
  edm::InputTag recoJets_;
  // output level manager
  int verbose;
  // event weight
  edm::InputTag weight_;
  // create gen plots?
  bool genPlots_;
  // within visible cross section range?
  edm::InputTag ingenPS_;
  // create rec plots?
  bool recPlots_;
  // use b-jet kinematics as before the kinematic fit?
  bool preBjets;
  // from gen level b-jet identification  
  edm::InputTag bHadJetIdx_, antibHadJetIdx_;
  // take minimum dR as criteria to assign rec and gen level b-jets to each other?
  // this is useful for simultaneous measurements of b-jet and anti-bjet 
  // where only the resolution and not the association matters
  bool useClosestDrBs_;

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
  float valueBqPtRec;
  float valueBqPtGen;
  float valueLeadBqPtRec;
  float valueLeadBqPtGen;
  float valueSubLeadBqPtRec;
  float valueSubLeadBqPtGen;
  float valueBqEtaRec;
  float valueBqEtaGen;
  float valueBqYRec;
  float valueBqYGen;
  float valueBbarqPtRec;
  float valueBbarqPtGen;
  float valueBbarqEtaRec;
  float valueBbarqEtaGen;
  float valueBbarqYRec;
  float valueBbarqYGen;
  float valueBbbarPtRec;
  float valueBbbarPtGen;
  float valueBbbarEtaRec;
  float valueBbbarEtaGen;
  float valueBbbarYRec;
  float valueBbbarYGen;
  float valueBbbarMassRec;
  float valueBbbarMassGen;
  float valueLbMassRec;
  float valueLbMassGen;
  float valueLepPtRec;
  float valueLepPtGen;
  float valueLepEtaRec;
  float valueLepEtaGen;
  float valueLepYRec;
  float valueLepYGen;

  //float valuexBHadRec;
  //float valuexBHadGen;
  //float valuexBLepRec;
  //float valuexBLepGen;
  bool bbSwapBetter;
  float valueAssignment;
  bool inVisPS;

  // ---
  //    histos
  // ---
  TH1F* bqPtRec;
  TH1F* bqPtGen;

  TH1F* bqPtLeadRec;
  TH1F* bqPtLeadGen;

  TH1F* bqPtSubLeadRec;
  TH1F* bqPtSubLeadGen;

  TH1F* bqEtaRec;
  TH1F* bqEtaGen;

  TH1F* bqYRec;
  TH1F* bqYGen;

  TH1F* bbbarPtRec;
  TH1F* bbbarPtGen;

  TH1F* bbbarEtaRec;
  TH1F* bbbarEtaGen;

  TH1F* bbbarYRec;
  TH1F* bbbarYGen;

  TH1F* bbbarMassRec;
  TH1F* bbbarMassGen;

  TH1F* lbMassRec;
  TH1F* lbMassGen;

  //TH1F* xBLepRec;
  //TH1F* xBLepGen;
  //TH1F* xBHadRec;
  //TH1F* xBHadGen;
  //TH1F* xBRec;
  //TH1F* xBGen;

  TH2F* bqPt_;
  TH2F* bqPtLead_;
  TH2F* bqPtSubLead_;
  TH2F* bqEta_;
  TH2F* bqY_;
  TH2F* bqPtClosestPt_;
  TH2F* bqEtaClosestPt_;
  TH2F* bqYClosestPt_;
  TH2F* bbbarPt_;
  TH2F* bbbarEta_;
  TH2F* bbbarY_;
  TH2F* bbbarMass_;
  TH2F* lbMass_;
  //TH2F* xB_;
};

#endif
