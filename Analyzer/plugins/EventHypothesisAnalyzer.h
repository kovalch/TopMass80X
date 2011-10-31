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
	edm::InputTag noPtEtaJets_;
  edm::InputTag leps_;
  edm::InputTag VertexSrc_;
  edm::InputTag PUWeightSrc_;
	edm::InputTag PUWeightUpSrc_;
	edm::InputTag PUWeightDownSrc_;
  edm::InputTag bWeightSrc_;
  edm::InputTag bWeightSrc_bTagSFUp_;
  edm::InputTag bWeightSrc_bTagSFDown_;
  edm::InputTag bWeightSrc_misTagSFUp_;
  edm::InputTag bWeightSrc_misTagSFDown_;
  edm::InputTag muWeightSrc_;
  
  bool savePDFWeights_;

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
  
  double hadQRawPt;
  double genHadQPt;
  
  double hadQBarPt;
  double hadQBarEta;
  double hadQBarMass;
  double hadQBarE;
  double hadQBarBSSV;
  double hadQBarJC;
  
  double hadQBarRawPt;
  
  double hadWPt;
  double hadWEta;
  double hadWMass;
  double hadWE;
  
  double hadWRawMass;
  double hadWRawPt;
  
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
  
  double hadBRawPt;
  
  double leptonPt;
  double leptonC;
	double leptonEta;
  
  double leptonRawPt;
  
  double neutrinoPt;
  double neutrinoE;
  
  double neutrinoRawPt;
  
  double lepWPt;
  double lepWEta;
  double lepWMass;
  double lepWE;
  
  double lepWRawMass;
  
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
  double hadTopRawPt;
  double hadTopRawEta;
  
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
  double deltaRHadBLepB;
  double deltaThetaHadBLepB;
  
  ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double> > sumB;
  double sumBPt;
  
  ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double> > TTBar;
  double TTBarPt;
  double TTBarMass;
  
  int jetMultiplicity;
  int noPtEtaJetMultiplicity;
  
  double noPtEtaJetPt;
  double leadingJetPt;
  
  double nlJetPt;
  double nlJetEta;
  
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
  double cProb;
  
  int nVertex;
  double PUWeight;
	double PUWeightUp;
	double PUWeightDown;
  double bWeight;
  double bWeight_bTagSFUp;
  double bWeight_bTagSFDown;
  double bWeight_misTagSFUp;
  double bWeight_misTagSFDown;
  double muWeight;
  double MCWeight;
  double pdfWeights[44];
  
  int target;
  
  TTree* eventTree;
  
  double QBTagProbabilitySSV(double bDiscriminator);
  double QBTagProbabilitySSVHEM(double bDiscriminator);
};

#endif
