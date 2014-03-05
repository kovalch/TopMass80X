#ifndef MixedObjectsAnalyzer_h
#define MixedObjectsAnalyzer_h

#include <map>
#include <vector>
#include <string>

#include "TTree.h"
#include "TH1F.h"
#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"
#include "AnalysisDataFormats/TopObjects/interface/TtSemiLepEvtPartons.h"
#include "AnalysisDataFormats/TopObjects/interface/TtSemiLeptonicEvent.h"

/**
   \class   MixedObjectsAnalyzer MixedObjectsAnalyzer.h "TopAnalysis/TopAnalyzer/interface/MixedObjectsAnalyzer.h"

   \brief   EDAnalyzer dealing with variables computed from different object types

   Computes quantities arising from different analysis objects such as invariant masses of jets, leptons and MET, etc. save them within a tree as well as histogram for flexible analysis input
*/

class MixedObjectsAnalyzer : public edm::EDAnalyzer {

 public:
  /// default constructor
  explicit MixedObjectsAnalyzer(const edm::ParameterSet& cfg);
  /// default destructor
  ~MixedObjectsAnalyzer();

 private:

  /// initiate histograms
  virtual void beginJob();
  /// produce n-tuple
  virtual void analyze(const edm::Event& event, const edm::EventSetup& iSetup);
  /// empty
  virtual void endJob();
  /// helper function to deal with abs(deltaphi) > pi
  double modTwoPi(double DeltaPhi);
  /// helper function to find to find non-ttbar jets
  bool findAncestor(const reco::Candidate* cand, TString decaychain);
  // helper function to get the jet-parton matching result for the KinFit reco jet assignment
  double checkPartonAssignment(const edm::Handle<TtSemiLeptonicEvent> semiLepEvt, int maxNJets);

  /// src's for the different infos
  edm::InputTag JetSrc_, METSrc_, MuonSrc_, ElectronSrc_,  GenJetSrc_, addGenJetSrc_, GenMETSrc_, GenLepSrc_, ingenPS_, weight_, VertexSrc_, semiLepEvt_;

  // class key of kinfit hypothesis
  std::string hypoKey_, btagAlgo_;

  /// event weight
  float weight;

  /// event identifiers
  unsigned int runNumber, luminosityBlockNumber, eventNumber;

  /// define Tree for event content
  TTree * tree;

  /// floats
  float btagDiscr_, MuNu4J, ElNu4J, mJJ, mWJJ, mWFitJJ, mHbb, leadNonttjetPt, leadNonttjetY, leadNonttjetEta;
  float addJetPt_;
  float bqhadPtPre, bqhadEtaPre, bqhadPhiPre, bqlepPtPre, bqlepEtaPre, bqlepPhiPre, lqPtPre, lqEtaPre, lqPhiPre, lqbarPtPre, lqbarEtaPre, lqbarPhiPre, nuPtPre, nuEtaPre, nuPhiPre, lepPtPre, lepEtaPre, lepPhiPre;
  float bqhadPtFit, bqhadEtaFit, bqhadPhiFit, bqlepPtFit, bqlepEtaFit, bqlepPhiFit, lqPtFit, lqEtaFit, lqPhiFit, lqbarPtFit, lqbarEtaFit, lqbarPhiFit, nuPtFit, nuEtaFit, nuPhiFit, lepPtFit, lepEtaFit, lepPhiFit;
  float nuPtTrue, nuEtaTrue, nuPhiTrue, lepPtTrue, lepEtaTrue, lepPhiTrue, sumEtPre, sumEtTrue;
  float nPV;
  float ttbarJetMass, ttbarJetMassTrue, rhos, rhosTrue;
  float leadNonttjetPtTrue, leadNonttjetYTrue, leadNonttjetEtaTrue;
  float topPtLepFit, topPtHadFit, topYLepFit, topYHadFit, ttbarMassFit, ttbarYFit, ttbarPtFit;
  float topPtLepTrue, topPtHadTrue, topYLepTrue, topYHadTrue, ttbarMassTrue, ttbarYTrue, ttbarPtTrue;
  float Nbjets, Njets, NjetsTrue;
  bool inVisPS;
  float valueAssignment;
  float chi2;

  /// ints
  int BindexA, BindexB, BindexC, BindexD, leadNonttjetIX;

  /// PDF related
  float Q;
  int id1, id2;
  double x1, x2;

  /// histo container
  std::map< std::string, TH1F* > hists_;
  std::map< std::string, TH2F* > hists2D_;
};

#endif
