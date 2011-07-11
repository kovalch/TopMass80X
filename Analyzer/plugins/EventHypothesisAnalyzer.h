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
  edm::InputTag jets_;
  edm::InputTag leps_;

  int run;
  int luminosityBlock;
  int event;
  int combi;
  
  double hadQPt;
  double hadQEta;
  double hadQMass;
  double hadQE;
  double hadQBSSV;
  double hadQJC;
  
  double hadQRawE;
  
  double hadQBarPt;
  double hadQBarEta;
  double hadQBarMass;
  double hadQBarE;
  double hadQBarBSSV;
  double hadQBarJC;
  
  double hadQBarRawE;
  
  double hadWPt;
  double hadWEta;
  double hadWMass;
  double hadWE;
  
  double hadWRawMass;
  double hadWRawSigM;
  
  double genHadWPt;
  double genHadWEta;
  double genHadWMass;
  double genHadWE;
  
  double hadBPt;
  double hadBEta;
  double hadBMass;
  double hadBE;
  double hadBBSSV;
  double hadBJC;
  
  double hadBRawE;
  
  double leptonPt;
  double leptonC;
  
  double lepBPt;
  double lepBEta;
  double lepBMass;
  double lepBE;
  double lepBBSSV;
  double lepBJC;
  
  double lepBRawE;
  
  double genHadBPt;
  double genHadBEta;
  double genHadBMass;
  double genHadBE;
  
  double hadTopPt;
  double hadTopEta;
  double hadTopMass;
  double hadTopE;
  
  double hadTopRawMass;
  
  double lepTopPt;
  double lepTopEta;
  double lepTopMass;
  double lepTopE;
  
  double lepTopRawMass;
  
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
  double hitFitChi2;
  double hitFitProb;
  double hitFitMT;
  double hitFitSigMT;
  double bProb;
  double hadBProb;
  double bProbSSV;
  double hadBProbSSV;
  
  int target;
  
  TTree* eventTree;
  
  double QBTagProbabilitySSV(double bDiscriminator);
};

#endif
