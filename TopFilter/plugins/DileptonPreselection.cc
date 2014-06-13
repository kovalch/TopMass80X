// -*- C++ -*-
//
// Package:    DileptonPreselection
// Class:      DileptonPreselection
// 
/**\class DileptonPreselection DileptonPreselection.cc TopAnalysis/TopFilter/plugins/DileptonPreselection.cc

 Description: Simple preselection, requiring at least one dilepton pair fulfilling the specified criteria

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Johannes Hauk,01a/O2.115,3139,
//         Created:  Fri Apr 11 14:15:34 CEST 2014
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

//
// class declaration
//

class DileptonPreselection : public edm::EDFilter {
   public:
      explicit DileptonPreselection(const edm::ParameterSet&);
      ~DileptonPreselection();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      
      /// Enumeration for possible channel selections
      enum Channel{ee, emu, mumu, undefined};
      
      /// Convert the channel from type string to enum
      Channel setChannel(const std::string& channel);
      
      /// Check the validity of specified parameters
      void checkParameterValidity()const;
      
      /// Does the lepton pair survive the selection requirements
      bool isGoodPair(const reco::RecoCandidate* lepton1, const reco::RecoCandidate* lepton2)const;
      
      // ----------member data ---------------------------
      
      /// Electron collection input tag
      const edm::InputTag electrons_;
      
      /// Muon collection input tag
      const edm::InputTag muons_;
      
      /// Requirement for the charges of the leptons
      const int filterCharge_;
      
      /// Require a specific dilepton channel
      const Channel filterChannel_;
      
      /// Vector for defining intervals of dilepton masses to be excluded
      const std::vector<double> v_excludeMasses_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DileptonPreselection::DileptonPreselection(const edm::ParameterSet& iConfig):
electrons_(iConfig.getParameter<edm::InputTag>("electrons")),
muons_(iConfig.getParameter<edm::InputTag>("muons")),
filterCharge_(iConfig.getParameter<int>("filterCharge")),
filterChannel_(this->setChannel(iConfig.getParameter<std::string>("filterChannel"))),
v_excludeMasses_(iConfig.getParameter<std::vector<double> >("excludeMasses"))
{
    this->checkParameterValidity();
}



DileptonPreselection::~DileptonPreselection()
{}



//
// member functions
//

DileptonPreselection::Channel DileptonPreselection::setChannel(const std::string& channel)
{
    if(channel == "ee") return ee;
    if(channel == "emu") return emu;
    if(channel == "mumu") return mumu;
    if(channel == "") return undefined;
    edm::LogError("DileptonPreselection")<<"Error in config parameter 'filterChannel'. Invalid value: "<<channel
                                         <<"\n...break";
    throw cms::Exception("Configuration Error");
}



void DileptonPreselection::checkParameterValidity()const
{
    if(filterCharge_<-1 || filterCharge_>1){
        edm::LogError("DileptonPreselection")<<"Error in config parameter 'filterCharge'. Invalid value: "<<filterCharge_
                                             <<"\n...break";
        throw cms::Exception("Configuration Error");
    }
    
    if(v_excludeMasses_.size()%2 == 1){
        edm::LogError("DileptonPreselection")<<"Error in config parameter 'excludeMasses'. "
                                             <<"Number of arguments is odd, but needs to even for defining intervals";
        throw cms::Exception("Configuration Error");
    }
    
    int entry(1);
    double intervalBegin(9999.);
    for(std::vector<double>::const_iterator i_excludeMasses = v_excludeMasses_.begin(); i_excludeMasses != v_excludeMasses_.end(); ++i_excludeMasses, ++entry){
        if(entry%2 == 1) intervalBegin = *i_excludeMasses;
        if(entry%2==0 && intervalBegin>*i_excludeMasses){
            edm::LogError("DileptonPreselection")<<"Incorrect Interval Definition in config parameter 'excludeMasses': \t"
                                                 <<intervalBegin<<" is bigger than "<<*i_excludeMasses
                                                 <<" but is expected to be smaller\n...break";
            throw cms::Exception("Configuration Error");
        }
    }
}



bool DileptonPreselection::isGoodPair(const reco::RecoCandidate* lepton1, const reco::RecoCandidate* lepton2)const
{
    // Exclude pairs with wrong charge combination
    if(filterCharge_){
        const int chargeProduct = lepton1->charge()*lepton2->charge();
        if(chargeProduct != filterCharge_) return false;
    }
    
    // Exclude pairs with wrong dilepton invariant mass
    if(v_excludeMasses_.size()){
        int entry(1);
        double intervalBegin(9999.);
        const double dileptonMass = (lepton1->p4() + lepton2->p4()).M();
        for(std::vector<double>::const_iterator i_excludeMasses = v_excludeMasses_.begin(); i_excludeMasses != v_excludeMasses_.end(); ++i_excludeMasses, ++entry){
            if(entry%2 == 1) intervalBegin = *i_excludeMasses;
            if(entry%2 == 0){
                if(dileptonMass>=intervalBegin && dileptonMass<*i_excludeMasses) return false;
            }
        }
    }
    
    // The pair survived, so a good lepton pair exits
    return true;
}



// ------------ method called on each new Event  ------------
bool
DileptonPreselection::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    // Access electrons and muons
    edm::Handle<std::vector<pat::Electron> > electrons;
    iEvent.getByLabel(electrons_, electrons);
    edm::Handle<std::vector<pat::Muon> > muons;
    iEvent.getByLabel(muons_, muons);
    
    size_t nElectrons = electrons->size();
    size_t nMuons = muons->size();
    
    // Skip events with less than 2 relevant leptons
    if(filterChannel_==undefined && nElectrons+nMuons<2) return false;
    if(filterChannel_==ee && nElectrons<2) return false;
    if(filterChannel_==emu && !nElectrons && !nMuons) return false;
    if(filterChannel_==mumu && nMuons<2) return false;
    
    // Check all ee candidates
    if(filterChannel_==undefined || filterChannel_==ee){
        for(std::vector<pat::Electron>::const_iterator i_electron = electrons->begin(); i_electron != electrons->end(); ++i_electron){
            for(std::vector<pat::Electron>::const_iterator j_electron = i_electron+1; j_electron != electrons->end(); ++j_electron){
                if(this->isGoodPair(&*i_electron, &*j_electron)) return true;
            }
        }
    }
    
    // Check all emu and mumu candidates
    for(std::vector<pat::Muon>::const_iterator i_muon = muons->begin(); i_muon != muons->end(); ++i_muon){
        // Check all emu candidates
        if(filterChannel_==undefined || filterChannel_==emu){
            for(std::vector<pat::Electron>::const_iterator i_electron = electrons->begin(); i_electron != electrons->end(); ++i_electron){
                if(this->isGoodPair(&*i_muon, &*i_electron)) return true;
            }
        }
        // Check all mumu candidates
        if(filterChannel_==undefined || filterChannel_==mumu){
            for(std::vector<pat::Muon>::const_iterator j_muon = i_muon+1; j_muon != muons->end(); ++j_muon){
                if(this->isGoodPair(&*i_muon, &*j_muon)) return true;
            }
        }
    }
    
    // No good pair is found
    return false;
}



// ------------ method called once each job just before starting event loop  ------------
void 
DileptonPreselection::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DileptonPreselection::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
DileptonPreselection::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
DileptonPreselection::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
DileptonPreselection::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
DileptonPreselection::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DileptonPreselection::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    // Prefer to set the parameters without a default value:
    // An exception is thrown when the parameter is not defined in the config files, instead of silently using the default given here
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("electrons");
    desc.add<edm::InputTag>("muons");
    desc.add<int>("filterCharge");
    desc.add<std::string>("filterChannel");
    desc.add<std::vector<double> >("excludeMasses");
    descriptions.add("dileptonPreselection", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DileptonPreselection);
