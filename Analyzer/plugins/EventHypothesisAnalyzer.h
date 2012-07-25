#ifndef EventHypothesisAnalyzer_h
#define EventHypothesisAnalyzer_h

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

#include "Math/VectorUtil.h"
#include "Math/GenVector/LorentzVector.h"
#include "TLorentzVector.h"
#include <Math/Boost.h>

#include "FWCore/Framework/interface/EDAnalyzer.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >	PxPyPzEVector;
typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >	PtEtaPhiMVector;

class EventHypothesisAnalyzer : public edm::EDAnalyzer {

 public:

  explicit EventHypothesisAnalyzer(const edm::ParameterSet&);
  ~EventHypothesisAnalyzer();
  
 private:

  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  edm::InputTag semiLepEvt_;
  edm::InputTag hypoClassKey_;
  
  edm::InputTag jets_;
  edm::InputTag allJets_;
  edm::InputTag noPtEtaJets_;
  edm::InputTag leps_;
  edm::InputTag mets_;

  edm::InputTag PUSrc_;
  edm::InputTag VertexSrc_;

  edm::InputTag PUWeightSrc_;
  edm::InputTag PUWeightUpSrc_;
  edm::InputTag PUWeightDownSrc_;

  edm::InputTag PUAWeightSrc_;
  edm::InputTag PUAWeightUpSrc_;
  edm::InputTag PUAWeightDownSrc_;

  edm::InputTag PUBWeightSrc_;
  edm::InputTag PUBWeightUpSrc_;
  edm::InputTag PUBWeightDownSrc_;

  edm::InputTag PUABWeightSrc_;
  edm::InputTag PUABWeightUpSrc_;
  edm::InputTag PUABWeightDownSrc_;

  edm::InputTag bWeightSrc_;
  edm::InputTag bWeightSrc_bTagSFUp_;
  edm::InputTag bWeightSrc_bTagSFDown_;
  edm::InputTag bWeightSrc_misTagSFUp_;
  edm::InputTag bWeightSrc_misTagSFDown_;
  
  edm::InputTag muWeightSrc_;
  edm::InputTag mcWeightSrc_;
  
  bool savePDFWeights_;
  bool data_;

  int run;
  int luminosityBlock;
  int event;
  int combi;
  
  double hadQPt;
  double hadQEta;
  double hadQMass;
  double hadQE;
  double hadQBSSV;
  double hadQBCSV;
  double hadQJC;
  double hadQF;
  
  double hadQRawPt;
  double hadQGenPt;
  
  double hadQBarPt;
  double hadQBarEta;
  double hadQBarMass;
  double hadQBarE;
  double hadQBarBSSV;
  double hadQBarBCSV;
  double hadQBarJC;
  double hadQBarF;
  
  double hadQBarRawPt;
  
  double hadWPt;
  double hadWEta;
  double hadWMass;
  double hadWE;
  
  double hadWRawMass;
  double hadWRawMass0;
  double hadWRawPt;
  
  double hadWGenPt;
  double hadWGenEta;
  double hadWGenMass;
  double hadWGenE;
  
  double hadBPt;
  double hadBEta;
  double hadBMass;
  double hadBE;
  double hadBBSSV;
  double hadBBCSV;
  double hadBJC;
  double hadBF;
  
  double hadBRawPt;
  
  double leptonPt;
  double leptonC;
  double leptonEta;
  int leptonId;
  
  double leptonRawPt;
  
  double nuPt;
  double nuE;
  
  double nuRawPt;
  
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
  double lepBBCSV;
  double lepBJC;
  double lepBF;
  
  double lepBRawE;
  
  double hadBGenPt;
  double hadBGenEta;
  double hadBGenMass;
  double hadBGenE;
  
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
  
  double hadTopGenPt;
  double hadTopGenEta;
  double hadTopGenMass;
  
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
  double jetsBCSV[6];
  int permutation[4];
  int noPtEtaJetMultiplicity;
  int bottomSSVJetMultiplicity;
  int bottomCSVJetMultiplicity;
  double jetsPt [6];
  double jetsEta[6];
  double jetsPhi[6];

  int allJetMultiplicity;
  int inWJetMultiplicity;
  double inWJetsPt [100];
  double inWJetsEta[100];
  double inWJetsPhi[100];
  double inWJetsDRQ[100];
  double inWJetsDRQBar[100];
  double inWJetsCHM[100];
  double inWJetsNC [100];
  int inWContamination;
  
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
  
  double amwtHadQ;
  double amwtHadQBar;
  double amwtLepton;
  
  double inWEta;
  double inWPhi;
  double inWdR;
  
  int nPU[3];
  int nVertex;
  
  double PUWeight;
  double PUWeightUp;
  double PUWeightDown;

  double PUAWeight;
  double PUAWeightUp;
  double PUAWeightDown;

  double PUBWeight;
  double PUBWeightUp;
  double PUBWeightDown;

  double PUABWeight;
  double PUABWeightUp;
  double PUABWeightDown;

  double bWeight;
  double bWeight_bTagSFUp;
  double bWeight_bTagSFDown;
  double bWeight_misTagSFUp;
  double bWeight_misTagSFDown;
  
  double muWeight;
  double mcWeight;
  
  double MCWeight;
  
  double pdfWeights[44];
  
  int target;
  
  TTree* eventTree;
  
  double QBTagProbabilitySSV(double bDiscriminator);
  double QBTagProbabilitySSVHEM(double bDiscriminator);
  
  PxPyPzEVector hadTop_;
  PxPyPzEVector hadTopRaw_;
  PxPyPzEVector hadTopGen_;
  
  PxPyPzEVector hadB_;
  PxPyPzEVector hadBRaw_;
  PxPyPzEVector hadBGen_;
  
  PxPyPzEVector hadW_;
  PxPyPzEVector hadWRaw_;
  PxPyPzEVector hadWGen_;
  PtEtaPhiMVector hadWGeom_;
  
  PxPyPzEVector hadQ_;
  PxPyPzEVector hadQRaw_;
  PxPyPzEVector hadQGen_;
  
  PxPyPzEVector hadQBar_;
  PxPyPzEVector hadQBarRaw_;
  PxPyPzEVector hadQBarGen_;
  
  PxPyPzEVector lepTop_;
  PxPyPzEVector lepTopRaw_;
  PxPyPzEVector lepTopGen_;
  
  PxPyPzEVector lepB_;
  PxPyPzEVector lepBRaw_;
  PxPyPzEVector lepBGen_;
  
  PxPyPzEVector lepW_;
  PxPyPzEVector lepWRaw_;
  PxPyPzEVector lepWGen_;
  
  PxPyPzEVector lepton_;
  PxPyPzEVector leptonRaw_;
  PxPyPzEVector leptonGen_;
  
  PxPyPzEVector nu_;
  PxPyPzEVector nuRaw_;
  PxPyPzEVector nuGen_;
};

#endif
