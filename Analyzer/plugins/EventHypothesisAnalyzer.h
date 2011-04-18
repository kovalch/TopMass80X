#ifndef EventHypothesisAnalyzer_h
#define EventHypothesisAnalyzer_h

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

#include "Math/VectorUtil.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"

class EventHypothesisAnalyzer : public edm::EDAnalyzer {

 public:

  explicit EventHypothesisAnalyzer(const edm::ParameterSet&);
  ~EventHypothesisAnalyzer(){};
  
 private:

  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  edm::InputTag semiLepEvt_;
  edm::InputTag hypoClassKey_;

  int run;
  int luminosityBlock;
  int event;
  int combi;
  
  double hadQPt;
  double hadQEta;
  double hadQMass;
  double hadQE;
  double hadQBTCHE;
  double hadQBVMVA;
  
  double hadQbarPt;
  double hadQbarEta;
  double hadQbarMass;
  double hadQbarE;
  double hadQbarBTCHE;
  double hadQbarBVMVA;
  
  double hadWPt;
  double hadWEta;
  double hadWMass;
  double hadWE;
  
  double genHadWPt;
  double genHadWEta;
  double genHadWMass;
  double genHadWE;
  
  double hadBPt;
  double hadBEta;
  double hadBMass;
  double hadBE;
  double hadBBTCHE;
  double hadBBVMVA;
  
  double lepBPt;
  double lepBEta;
  double lepBMass;
  double lepBE;
  double lepBBTCHE;
  double lepBBVMVA;
  
  double genHadBPt;
  double genHadBEta;
  double genHadBMass;
  double genHadBE;
  
  double hadTopPt;
  double hadTopEta;
  double hadTopMass;
  double hadTopE;
  
  double genHadTopPt;
  double genHadTopEta;
  double genHadTopMass;
  
  double deltaRHadQHadQBar;
  double deltaThetaHadQHadQBar;
  double deltaRHadWHadB;
  double deltaThetaHadWHadB;
  double genDeltaThetaHadWHadB;
  double deltaRLepBLepton;
  double deltaThetaLepBLepton;
  
  double genMatchDr;
  double mvaDisc;
  double fitChi2;
  double fitProb;
  
  int target;
  
  TTree* eventTree;

};

#endif
