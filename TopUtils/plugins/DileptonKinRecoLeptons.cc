// -*- C++ -*-
//
// Package:    DileptonKinRecoLeptons
// Class:      DileptonKinRecoLeptons
// 
/**\class DileptonKinRecoLeptons DileptonKinRecoLeptons.cc TopAnalysis/TopUtils/plugins/DileptonKinRecoLeptons.cc

 Description: Returns the two leptons which should be fed in the dilepton kinematic reconstruction

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Johannes Hauk,01a/O2.115,3139,
//         Created:  Sat Apr 12 14:12:06 CEST 2014
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"


//
// class declaration
//

class DileptonKinRecoLeptons : public edm::EDProducer {
   public:
      explicit DileptonKinRecoLeptons(const edm::ParameterSet&);
      ~DileptonKinRecoLeptons();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      
      /// Enumeration for possible channel selections
      enum Channel{ee, emu, mumu, undefined};
      
      /// Convert the channel from type string to enum
      Channel setChannel(const std::string& channel);
      
      /// Check the validity of specified parameters
      void checkParameterValidity()const;
      
      /// Does the lepton pair survive the selection requirements
      bool isGoodPair(const reco::RecoCandidate* lepton1, const reco::RecoCandidate* lepton2)const;
      
      /// Orders the lepton pair by pt
      std::pair<const reco::RecoCandidate*, const reco::RecoCandidate*> ptOrderedPair(const reco::RecoCandidate* lepton1, const reco::RecoCandidate* lepton2)const;
      
      // ----------member data ---------------------------
      
      /// Electron collection input tag
      const edm::InputTag electrons_;
      
      /// Muon collection input tag
      const edm::InputTag muons_;
      
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
DileptonKinRecoLeptons::DileptonKinRecoLeptons(const edm::ParameterSet& iConfig):
electrons_(iConfig.getParameter<edm::InputTag>("electrons")),
muons_(iConfig.getParameter<edm::InputTag>("muons")),
filterChannel_(this->setChannel(iConfig.getParameter<std::string>("filterChannel"))),
v_excludeMasses_(iConfig.getParameter<std::vector<double> >("excludeMasses"))
{
    produces<std::vector<pat::Electron> >();
    produces<std::vector<pat::Muon> >();
}


DileptonKinRecoLeptons::~DileptonKinRecoLeptons()
{}


//
// member functions
//

DileptonKinRecoLeptons::Channel DileptonKinRecoLeptons::setChannel(const std::string& channel)
{
    if(channel == "ee") return ee;
    if(channel == "emu") return emu;
    if(channel == "mumu") return mumu;
    if(channel == "") return undefined;
    edm::LogError("DileptonKinRecoLeptons")<<"Error in config parameter 'filterChannel'. Invalid value: "<<channel
                                         <<"\n...break";
    throw cms::Exception("Configuration Error");
}



void DileptonKinRecoLeptons::checkParameterValidity()const
{
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



bool DileptonKinRecoLeptons::isGoodPair(const reco::RecoCandidate* lepton1, const reco::RecoCandidate* lepton2)const
{
    // Exclude pairs with wrong charge combination
    if(lepton1->charge() == lepton2->charge()) return false;
    
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



std::pair<const reco::RecoCandidate*, const reco::RecoCandidate*> DileptonKinRecoLeptons::ptOrderedPair(const reco::RecoCandidate* lepton1,
                                                                                                        const reco::RecoCandidate* lepton2)const
{
    if(!lepton1 || !lepton2){
        edm::LogError("DileptonPreselection")<<"Error in ptOrderedPair(). At least one of the leptons is a null pointer\n...break";
        throw cms::Exception("Wrong input for method");
    }
    
    const reco::RecoCandidate* leadingLepton(0);
    const reco::RecoCandidate* subleadingLepton(0);
    if(lepton1->pt() >= lepton2->pt()){
        leadingLepton = &*lepton1;
        subleadingLepton = &*lepton2;
    }
    else{
        leadingLepton = &*lepton2;
        subleadingLepton = &*lepton1;
    }
    
    return std::make_pair(leadingLepton, subleadingLepton);
}



// ------------ method called to produce the data  ------------
void
DileptonKinRecoLeptons::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    // Result vectors for the two electrons+muons
    std::auto_ptr<std::vector<pat::Electron> > selectedElectrons(new std::vector<pat::Electron>);
    std::auto_ptr<std::vector<pat::Muon> > selectedMuons(new std::vector<pat::Muon>);
    
    // Access electrons and muons
    edm::Handle<std::vector<pat::Electron> > electrons;
    iEvent.getByLabel(electrons_, electrons);
    size_t nElectrons = electrons->size();
    edm::Handle<std::vector<pat::Muon> > muons;
    iEvent.getByLabel(muons_, muons);
    size_t nMuons = muons->size();
    
    // Skip events with less than 2 relevant leptons, but create empty collections
    if((filterChannel_==undefined && nElectrons+nMuons<2) || (filterChannel_==ee && nElectrons<2) ||
       (filterChannel_==emu && !nElectrons && !nMuons) || (filterChannel_==mumu && nMuons<2)){
        iEvent.put(selectedElectrons);
        iEvent.put(selectedMuons);
        return;
    }
    
    // Vector containing all lepton pairs fulfilling the selection
    std::auto_ptr<std::vector<std::pair<const reco::RecoCandidate*, const reco::RecoCandidate*> > > leptonPairs(new std::vector<std::pair<const reco::RecoCandidate*, const reco::RecoCandidate*> >);
    
    // Check all ee candidates
    if(filterChannel_==undefined || filterChannel_==ee){
        for(std::vector<pat::Electron>::const_iterator i_electron = electrons->begin(); i_electron != electrons->end(); ++i_electron){
            for(std::vector<pat::Electron>::const_iterator j_electron = i_electron+1; j_electron != electrons->end(); ++j_electron){
                if(this->isGoodPair(&*i_electron, &*j_electron)) leptonPairs->push_back(this->ptOrderedPair(&*i_electron, &*j_electron));
            }
        }
    }
    
    // Check all emu and mumu candidates
    for(std::vector<pat::Muon>::const_iterator i_muon = muons->begin(); i_muon != muons->end(); ++i_muon){
        // Check all emu candidates
        if(filterChannel_==undefined || filterChannel_==emu){
            for(std::vector<pat::Electron>::const_iterator i_electron = electrons->begin(); i_electron != electrons->end(); ++i_electron){
                if(this->isGoodPair(&*i_muon, &*i_electron)) leptonPairs->push_back(this->ptOrderedPair(&*i_muon, &*i_electron));
            }
        }
        // Check all mumu candidates
        if(filterChannel_==undefined || filterChannel_==mumu){
            for(std::vector<pat::Muon>::const_iterator j_muon = i_muon+1; j_muon != muons->end(); ++j_muon){
                if(this->isGoodPair(&*i_muon, &*j_muon)) leptonPairs->push_back(this->ptOrderedPair(&*i_muon, &*j_muon));
            }
        }
    }
    
    // Find best pair as: search leading lepton, in case of several combinations choose the 2nd lepton with highest pt
    const reco::RecoCandidate* leadingLepton(0);
    const reco::RecoCandidate* subleadingLepton(0);
    for(std::vector<std::pair<const reco::RecoCandidate*, const reco::RecoCandidate*> >::const_iterator leptonPair = leptonPairs->begin(); leptonPair != leptonPairs->end(); ++leptonPair){
        const reco::RecoCandidate* lepton1(leptonPair->first);
        const reco::RecoCandidate* lepton2(leptonPair->second);
        if(!leadingLepton){
            leadingLepton = lepton1;
            subleadingLepton = lepton2;
        }
        else{
            if(leadingLepton == lepton1){
                if(lepton2->pt() > subleadingLepton->pt()) subleadingLepton = lepton2;
            }
            else if(lepton1->pt() > leadingLepton->pt()){
                leadingLepton = lepton1;
                subleadingLepton = lepton2;
            }
        }
    }
    
    // Fill result vectors for the two electrons+muons (from all dileptons the combination of: leading lepton, and corresponding most subleading lepton)
    if(leadingLepton && subleadingLepton){
        const pat::Electron* electron(0);
        electron = dynamic_cast<const pat::Electron*>(leadingLepton);
        if(electron) selectedElectrons->push_back(*electron);
        else selectedMuons->push_back(*static_cast<const pat::Muon*>(leadingLepton));
        electron = dynamic_cast<const pat::Electron*>(subleadingLepton);
        if(electron) selectedElectrons->push_back(*electron);
        else selectedMuons->push_back(*static_cast<const pat::Muon*>(subleadingLepton));
    }
    iEvent.put(selectedElectrons);
    iEvent.put(selectedMuons);
}

// ------------ method called once each job just before starting event loop  ------------
void 
DileptonKinRecoLeptons::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DileptonKinRecoLeptons::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
DileptonKinRecoLeptons::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
DileptonKinRecoLeptons::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
DileptonKinRecoLeptons::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
DileptonKinRecoLeptons::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DileptonKinRecoLeptons::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    // Prefer to set the parameters without a default value:
    // An exception is thrown when the parameter is not defined in the config files, instead of silently using the default given here
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("electrons");
    desc.add<edm::InputTag>("muons");
    desc.add<std::string>("filterChannel");
    desc.add<std::vector<double> >("excludeMasses");
    descriptions.add("dileptonKinRecoLeptons", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DileptonKinRecoLeptons);
