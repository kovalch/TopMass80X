#ifndef BRegJetEventAnalyzer_h
#define BRegJetEventAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "TopMass/TopEventTree/interface/BRegJetEvent.h"

// for event setup handles needed for "computer"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"


//for "computer" that is needed for more advanced b-tagging variables
#include "RecoBTau/JetTagComputer/interface/GenericMVAJetTagComputer.h"
#include "RecoBTau/JetTagComputer/interface/GenericMVAJetTagComputerWrapper.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputer.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputerRecord.h"
#include "RecoBTag/SecondaryVertex/interface/CombinedSVComputer.h"

//to read in regression results
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

// PUJetIdIdentifier
#include "CMGTools/External/interface/PileupJetIdentifier.h"
#include "CMGTools/External/interface/PileupJetIdAlgo.h"

//GBRTrain
#include "CondFormats/EgammaObjects/interface/GBRForest.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TopMass/TopEventTree/interface/TreeRegistryService.h"

#include "TFile.h"

//#include <utility>

class BRegJetEventAnalyzer : public edm::EDAnalyzer {

  public:

  explicit BRegJetEventAnalyzer(const edm::ParameterSet&);
  ~BRegJetEventAnalyzer();

  void fillBRegJet(const edm::Event&, const edm::EventSetup&);
  std::vector<double> fillBRegJetAndReturnGBRResults(const edm::Event&, const edm::EventSetup&);
  std::vector<double> fillBRegJetAndReturnTMVAResults(const edm::Event&, const edm::EventSetup&);
  void iniBRegEvent();

 private:

  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  edm::Service<TreeRegistryService> trs;

  edm::InputTag jets_;
  //edm::InputTag allJets_;
  //edm::InputTag noPtEtaJets_;
  edm::InputTag rho_tag_;
  edm::InputTag rho25_tag_;
  std::string gluonTagName_;

  // max possible number of jets in events
  const int kJetMAX_;

  // THE BRegJetEvent to store the information
  BRegJetEvent* BRegJet;

  //TMVA reader to determine b-regression result
  TMVA::Reader *reader_;
  std::string mva_name_,mva_path_,GBRmva_path_;



  bool writeOutVariables_;

  /// JECUncertaintySource File
  edm::FileInPath JECUncSrcFile_;

  /// jet output collection 
  std::string outputJets_;

  // check only once per module run if the needed collections are available
  bool checkedIsPFJet, checkedJERSF, checkedJESSF, checkedTotalSF, checkedQGTag, checkedJESTotUnc,  checkedBRegResult;
  bool        isPFJet,     hasJERSF,     hasJESSF,     hasTotalSF,     hasQGTag,     hasJESTotUnc,      hasBRegResult;

  double tempJetPtCorr_, tempJetEta_;

  std::map<std::string, double> AllVarHolder_;
  std::vector <std::string> AllVariables_;

  GBRForest *gbropt_;
  std::vector <double*> variablePointer_;
  Float_t *vals_;
  std::vector<std::string> *varlist_;

};

#endif
