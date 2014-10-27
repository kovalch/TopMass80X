// -*- C++ -*-
//
// Package:    LeptonVertexSelector
// Class:      LeptonVertexSelector
// 
/**\class LeptonVertexSelector LeptonVertexSelector.cc TopAnalysis/TopUtils/plugins/LeptonVertexSelector.cc

 Description: Select electrons and muons having maximum distances in dxy and dz with respect to the leading primary vertex

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Johannes Hauk,01a/O2.115,3139,
//         Created:  Fri Jul 18 13:55:44 CEST 2014
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

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"



//
// class declaration
//

class LeptonVertexSelector : public edm::EDProducer{
public:
    explicit LeptonVertexSelector(const edm::ParameterSet&);
    ~LeptonVertexSelector();
    
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
private:
    virtual void beginJob();
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob();
    
    virtual void beginRun(edm::Run&, edm::EventSetup const&);
    virtual void endRun(edm::Run&, edm::EventSetup const&);
    virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
    virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
    
    // ----------member data ---------------------------
    
    // Input collections
    const edm::InputTag electronsTag_;
    const edm::InputTag muonsTag_;
    const edm::InputTag verticesTag_;
    
    // Cut parameters
    const double electronDxyMax_;
    const double electronDzMax_;
    const double muonDxyMax_;
    const double muonDzMax_;
    
    // Should electrons/muons being cut?
    bool cutElectrons_;
    bool cutMuons_;
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
LeptonVertexSelector::LeptonVertexSelector(const edm::ParameterSet& iConfig):
electronsTag_(iConfig.getParameter<edm::InputTag>("electrons")),
muonsTag_(iConfig.getParameter<edm::InputTag>("muons")),
verticesTag_(iConfig.getParameter<edm::InputTag>("vertices")),
electronDxyMax_(iConfig.getParameter<double>("electronDxyMax")),
electronDzMax_(iConfig.getParameter<double>("electronDzMax")),
muonDxyMax_(iConfig.getParameter<double>("muonDxyMax")),
muonDzMax_(iConfig.getParameter<double>("muonDzMax")),
cutElectrons_(false),
cutMuons_(false)
{
    if(electronDxyMax_>=0. || electronDzMax_>=0.){
        cutElectrons_ = true;
        edm::LogInfo log("LeptonVertexSelector");
        log<<"Selection for electrons applied:\n";
        if(electronDxyMax_ >= 0.) log<<"\t |dxy| < "<<electronDxyMax_<<"\n";
        if(electronDzMax_ >= 0.) log<<"\t |dz| < "<<electronDzMax_<<"\n";
    }
    if(muonDxyMax_>=0. || muonDzMax_>=0.){
        cutMuons_ = true;
        edm::LogInfo log("LeptonVertexSelector");
        log<<"Selection for muons applied:\n";
        if(muonDxyMax_ >= 0.) log<<"\t |dxy| < "<<muonDxyMax_<<"\n";
        if(muonDzMax_ >= 0.) log<<"\t |dz| < "<<muonDzMax_<<"\n";
    }
    
    produces<std::vector<pat::Electron> >();
    produces<std::vector<pat::Muon> >();
}



LeptonVertexSelector::~LeptonVertexSelector()
{}



//
// member functions
//

// ------------ method called to produce the data  ------------
void
LeptonVertexSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    // Result vectors for the electrons and muons
    std::auto_ptr<std::vector<pat::Electron> > selectedElectrons(new std::vector<pat::Electron>);
    std::auto_ptr<std::vector<pat::Muon> > selectedMuons(new std::vector<pat::Muon>);
    
    // Access vertices
    edm::Handle<std::vector<reco::Vertex> > vertices;
    iEvent.getByLabel(verticesTag_, vertices);
    const size_t nVertex = vertices->size();
    if(!nVertex) edm::LogInfo("LeptonVertexSelector")<<"No primary vertex found in event";
    
    // Access electrons
    edm::Handle<std::vector<pat::Electron> > electrons;
    iEvent.getByLabel(electronsTag_, electrons);
    if(!cutElectrons_){
        selectedElectrons->assign(electrons->begin(), electrons->end());
    }
    else{
        if(nVertex){
            for(std::vector<pat::Electron>::const_iterator i_electron = electrons->begin(); i_electron != electrons->end(); ++i_electron){
                const reco::GsfTrack& track = *(i_electron->gsfTrack());
                if(electronDxyMax_ >= 0.){
                    const double dxyAbs = std::abs(track.dxy(vertices->at(0).position()));
                    if(dxyAbs > electronDxyMax_) continue;
                }
                if(electronDzMax_ >= 0.){
                    const double dzAbs = std::abs(track.dz(vertices->at(0).position()));
                    if(dzAbs > electronDzMax_) continue;
                }
                selectedElectrons->push_back(*i_electron);
            }
            //if(electrons->size()) std::cout<<"Electron sizes: "<<electrons->size()<<" , "<<selectedElectrons->size()<<"\n";
            //if(electrons->size() != selectedElectrons->size()) std::cout<<"Electrons: "<<electrons->size()<<" , "<<selectedElectrons->size()<<"\n";
        }
    }
    
    // Access muons
    edm::Handle<std::vector<pat::Muon> > muons;
    iEvent.getByLabel(muonsTag_, muons);
    if(!cutMuons_){
        selectedMuons->assign(muons->begin(), muons->end());
    }
    else{
        if(nVertex){
            for(std::vector<pat::Muon>::const_iterator i_muon = muons->begin(); i_muon != muons->end(); ++i_muon){
                if(muonDxyMax_ >= 0.){
                    const double dxyAbs = std::abs(i_muon->dB());
                    if(dxyAbs > muonDxyMax_) continue;
                }
                if(muonDzMax_ >= 0.){
                    reco::TrackRef track = i_muon->muonBestTrack();
                    if(!track.isAvailable()){
                        edm::LogWarning("LeptonVertexSelector")<<"Cannot find muon best track for dz cut!!!";
                        continue;
                    }
                    const double dzAbs = std::abs(track->dz(vertices->at(0).position()));
                    if(dzAbs > muonDzMax_) continue;
                }
                selectedMuons->push_back(*i_muon);
            }
            //if(muons->size()) std::cout<<"Muon sizes: "<<muons->size()<<" , "<<selectedMuons->size()<<"\n";
            //if(muons->size() != selectedMuons->size()) std::cout<<"Muons: "<<muons->size()<<" , "<<selectedMuons->size()<<"\n";
        }
    }
    
    iEvent.put(selectedElectrons);
    iEvent.put(selectedMuons);
}



// ------------ method called once each job just before starting event loop  ------------
void 
LeptonVertexSelector::beginJob()
{}



// ------------ method called once each job just after ending the event loop  ------------
void 
LeptonVertexSelector::endJob()
{}



// ------------ method called when starting to processes a run  ------------
void 
LeptonVertexSelector::beginRun(edm::Run&, edm::EventSetup const&)
{}



// ------------ method called when ending the processing of a run  ------------
void 
LeptonVertexSelector::endRun(edm::Run&, edm::EventSetup const&)
{}



// ------------ method called when starting to processes a luminosity block  ------------
void 
LeptonVertexSelector::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{}



// ------------ method called when ending the processing of a luminosity block  ------------
void 
LeptonVertexSelector::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{}



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LeptonVertexSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    // Prefer to set the parameters without a default value:
    // An exception is thrown when the parameter is not defined in the config files, instead of silently using the default given here
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("electrons");
    desc.add<edm::InputTag>("muons");
    desc.add<edm::InputTag>("vertices");
    desc.add<double>("electronDxyMax");
    desc.add<double>("electronDzMax");
    desc.add<double>("muonDxyMax");
    desc.add<double>("muonDzMax");
    descriptions.add("LeptonVertexSelector", desc);
}



//define this as a plug-in
DEFINE_FWK_MODULE(LeptonVertexSelector);
