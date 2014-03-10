#ifndef TopKinematics_h
#define TopKinematics_h

#include "DataFormats/Candidate/interface/Candidate.h"
#include "TopAnalysis/TopAnalyzer/interface/SingleObject.h"
#include "AnalysisDataFormats/TopObjects/interface/TtSemiLeptonicEvent.h"
#include "TMath.h"

/**
   \class   TopKinematics TopKinematics.h "TopAnalysis/TopAnalyzer/interface/TopKinematics.h"

   \brief   Derived class to analyze the kinematics of top quarks on reconstruction and generator level

   The structure keeps histograms for the kinematics of top quarks and ttbar pairs. These histograms 
   can be filled from the TtSemiLeptonicEvent class (with an extra option to require matching for 
   purity & stability studies) or from the TtGenEvent class. The class is derived from the 
   SingleObject<Collection> interface, which makes it usable in full framework or fwlite. 
*/
double pi = TMath::Pi();

class TopKinematics : public SingleObject<TtSemiLeptonicEvent> {
  
 public:
  /// default constructor for generator level analysis in fw lite
  explicit TopKinematics();
  /// default constructor for reco level analysis in fw lite
  explicit TopKinematics(const std::string& hypoKey, const std::string& lepton, const bool& matchForStabilityAndPurity);
  /// default constructor for fwfull
  explicit TopKinematics(const edm::ParameterSet& configFile);
  /// default destructor
  ~TopKinematics(){};

  /**
     The following functions have to be implemented for any class
     derived from SingleObject<Collection>. The fill function is
     overloaded to serve both for the TtGenEvent class directoy or 
     the TtSemiLeptonicEvent class.
  **/
  /// histogramm booking for fwlite 
  void book();
  /// histogramm booking for fw
  void book(edm::Service<TFileService>& fileService);
  /// histogram filling interface for reconstruction level for access with fw or fwlite
  void fill(const TtSemiLeptonicEvent& tops, const double& weight=1.);
  /// histogram filling interface for generator level for access with fw or fwlite 
  void fill(const TtGenEvent& tops, const double& weight=1.);
  // histogram filling and in addition save information for event identification
  // double definition to avoid problems with different module inputs
  void fill2(const TtGenEvent& tops, const double& runNumber, const double& luminosityBlockNumber, const double& eventNumber, const double& x1, const double& x2, const float& Q, const int& id1, const int& id2, const double& weight, std::vector<double> weights){
    fillValue("runNumber"            , runNumber            , weight);
    fillValue("luminosityBlockNumber", luminosityBlockNumber, weight);
    fillValue("eventNumber"          , eventNumber          , weight);
    fillValue("Q"  , Q  , weight);
    fillValue("id1", id1, weight);
    fillValue("id2", id2, weight);
    fillValue("x1" , x1 , weight);
    fillValue("x2" , x2 , weight);
     // if more than one weight is supposed to be stored
    for(unsigned int iWeight=0; iWeight < this->wgts_.size(); iWeight++){
      std::string weightName = this->wgts_[iWeight].label()+this->wgts_[iWeight].instance();
      if(this->wgts_.size() == this->weights.size()) {
	fillValue(weightName, this->weights[iWeight], weight);
      }
      else std::cout<< "ERROR!!! Size of weights ("<<this->weights.size()<<") != size of weight tags ("<<wgts_.size()<<")!!! No weights are filled in tree!" <<std::endl;
    }
    fill(tops, weight);
  }
  void fill2(const TtSemiLeptonicEvent& tops, const double& runNumber, const double& luminosityBlockNumber, const double& eventNumber, const double& x1, const double& x2, const float& Q, const int& id1, const int& id2, const double& weight, std::vector<double> weights){
    fillValue("runNumber", runNumber, weight);
    fillValue("luminosityBlockNumber", luminosityBlockNumber, weight);
    fillValue("eventNumber", eventNumber, weight);
    fillValue("id1", id1, weight);
    fillValue("id2", id2, weight);
    fillValue("x1" , x1 , weight);
    fillValue("x2" , x2 , weight);
    fillValue("Q"  , Q  , weight);
     // if more than one weight is supposed to be stored
    for(unsigned int iWeight=0; iWeight < this->wgts_.size(); iWeight++){
      std::string weightName = this->wgts_[iWeight].label()+this->wgts_[iWeight].instance();
      if(this->wgts_.size() == this->weights.size()) {
	fillValue(weightName, this->weights[iWeight], weight);
      }
      else std::cout<< "ERROR!!! Size of weights ("<<this->weights.size()<<") != size of weight tags ("<<wgts_.size()<<")!!! No weights are filled in tree!" <<std::endl;
    }
    fill(tops, weight);
  }
  /// everything which needs to be done after the event loop
  void process(){
/*     TFile &file = fs->file(); */
/*     for(std::map<std::string, TH2*>::const_iterator iter = corrs_.begin(); iter != corrs_.end(); ++iter){ */
/*       std::cout << self->moduleDescription_.moduleLabel() << std::endl; */
/*       std::cout << self->moduleDescription_.moduleName() << std::endl; */
/*       TString folder=""; */
/* 	//file.cd(folder); */
/*       if(!iter->second->GetEntries()){  */
/* 	std::cout << "deleting" << folder << "/" << iter->second->GetName() << std::endl; */
/* 	//file.Delete(TString(iter->second->GetName())+TString(";*")); */
/* 	std::cout << "ok" << std::endl; */
/*       } */
/*     } */
  };

 private:
  /**
     Helper functions for special histogram management
  **/
  /// histogram filling for candidates topA and topB
  void fill(const reco::Candidate* leptonicTop, const reco::Candidate* hadronicTop, 
	    const reco::Candidate* topPlus    , const reco::Candidate* topMinus   , 
	    const reco::Candidate* leptonicW  , const reco::Candidate* hadronicW  , 
	    const reco::Candidate* leptonicB  , const reco::Candidate* hadronicB  , 
	    double HT, const double& charge, const double& weight=1.);

  /// histogram filling for candidate topA and topB (for stability and purity calculation)
  void fillGenRec(const reco::Candidate* leptonicTopRec, const reco::Candidate* leptonicTopGen, 
	    const reco::Candidate* hadronicTopRec, const reco::Candidate* hadronicTopGen, 
	    const reco::Candidate* leptonicWRec  , const reco::Candidate* leptonicWGen  , 
	    const reco::Candidate* hadronicWRec  , const reco::Candidate* hadronicWGen  , 
	    double HTrec, double HTgen, const double& weight=1.);

  /// helper function for determining stability and purity
  void match(const std::string& histo, const double& genValue, const double& recValue, const double& weight);
  ///  helper functions to fill 1D angle histos
  void fillAngles(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& hadB    ,
		  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& q       ,
		  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& qbar    ,
		  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& lepB    ,
		  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& muon    ,
		  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& neutrino,
		  const double& weight=1.);
  ///  helper functions to fill event shape variables
  void fillEventShapes(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& hadB    ,
		       const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& q       ,
		       const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& qbar    ,
		       const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& lepB    ,
		       const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& muon    ,
		       const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& neutrino,
		       const double& weight, const bool useMu, const bool useNu);
  ///  helper functions to fill final state objects
  void fillFinalStateObjects(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& hadB    ,
			     const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& q       ,
			     const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& qbar    ,
			     const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& lepB    ,
			     const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& lepton  ,
			     const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& neutrino,
			     const double& charge,
			     const double& weight);
  /// helper function to get the decay mode
  double getDecayChannel(const TtGenEvent& tops);
  // identify event by the following number coding:
  // -4 : undefined
  // -3 : non ttbar MC
  // -2 : dileptonic unknown    ttbar MC
  // -1 : semileptonic unknown  ttbar MC
  // 01 : semileptonic electron ttbar MC
  // 02 : semileptonic muon     ttbar MC
  // 03 : semileptonic tau      ttbar MC
  // 11 : dileptonic e   e      ttbar MC
  // 12 : dileptonic e   mu     ttbar MC
  // 13 : dileptonic e   tau    ttbar MC
  // 14 : dileptonic mu  mu     ttbar MC
  // 15 : dileptonic mu  tau    ttbar MC
  // 16 : dileptonic tau tau    ttbar MC
  // 20 : fully hadronic        ttbar MC

  /// helper function to get the final state lepton in tau->lepton events 
  // decay chain is: tau(status 3)->tau(status 2)->e/mu(status 1)
  reco::GenParticle* getFinalStateLepton(const reco::GenParticle& particle);

 private:
  /// class key of hypothesis
  std::string hypoKey_;
  /// lepton flavour in semileptonic event: muon (default) or electron
  std::string lepton_;
  /// bool to decide whether to use a ttree 
  bool useTree_;
  /// apply matching for stability and purity or not
  bool matchForStabilityAndPurity_;
  /// choose whether you want to destinguish between top/antitop instead of leptonic/hadronic top
  bool ttbarInsteadOfLepHadTop_;
  /// number of jets considered in Kinematic fit (needed for the classification of permutations)
  int maxNJets;
  /// number degrees of freedom used to calculate fit probability
  int ndof;
  /// histogram container for correlation plots
  std::map<std::string, TH2*> corrs_;
};

#endif
