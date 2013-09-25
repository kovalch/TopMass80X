#ifndef DummyProducer_h
#define DummyProducer_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "TopMass/TopEventTree/interface/BRegJetEvent.h"
#include "TopMass/TopEventTree/plugins/BRegJetEventAnalyzer.h"

//for "computer" that is needed for more advanced b-tagging variables
#include "RecoBTau/JetTagComputer/interface/GenericMVAJetTagComputer.h"
#include "RecoBTau/JetTagComputer/interface/GenericMVAJetTagComputerWrapper.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputer.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputerRecord.h"
#include "RecoBTag/SecondaryVertex/interface/CombinedSVComputer.h"

//to read in regression results
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

//GBRTrain
#include "CondFormats/EgammaObjects/interface/GBRForest.h"


#include "TFile.h"

/**
   \class   DummyProducer DummyProducer.h "TopAnalysis/TopUtils/plugins/DummyProducer.h"

   \brief   Plugin to add b-jet specific data class to event and (if only jet-specific regression is selected)
   	   	   	scale jet energies

*/

class BRegProducer : public edm::EDProducer {

 public:
  /// default constructor
  explicit BRegProducer(const edm::ParameterSet&);
  /// default destructor
  ~BRegProducer(){};
  
 private:
  /// check settings
  virtual void beginJob();
  /// rescale jet energy and recalculated MET
  virtual void produce(edm::Event&, const edm::EventSetup&);
//  /// rescale the resolution of the jet
//  double resolutionFactor(const pat::Jet&);
//  /// scale all energies of the jet
//  void scaleJetEnergy(pat::Jet&, double);

 private:
  /// jet input collection 
  edm::InputTag inputJets_;
//  /// met input collection
//  edm::InputTag inputMETs_;

  /// jet output collection 
  std::string outputJets_;
//  /// MET output collection


  BRegJetEventAnalyzer *bRegAnalyzer;

};

#endif
