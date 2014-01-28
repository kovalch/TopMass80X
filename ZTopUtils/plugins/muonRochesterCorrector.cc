#include "FWCore/Utilities/interface/EDMException.h"
#include "muonRochesterCorrector.h"
#include "CommonTools/Utils/interface/PtComparator.h"

#include <vector>
#include <iostream>
#include <stdexcept>
#include <algorithm>

muonRochesterCorrector::muonRochesterCorrector(const edm::ParameterSet& cfg):
debug(cfg.getParameter< bool >("debug") ),
inputmuons_( cfg.getParameter<edm::InputTag>("muonSrc") ),
isMC_(cfg.getParameter< bool >("isMC") ),
muontypestring_(cfg.getParameter< std::string >("muonType")),
seed_(cfg.getParameter<  int >("randomSeed"))
{

    std::cout << "********* Initializing muonRochesterCorrector *********\n" ;
    std::cout << " -> Using muon type: " << muontypestring_;

    if(muontypestring_ == "patMuons"){
        produces<std::vector<pat::Muon> >();
        muontype_=patmuons;}
    else if (muontypestring_ == "recoMuons"){
        produces<std::vector<reco::Muon> >();
        muontype_=recomuons;}
    else if (muontypestring_ == "pfMuons"){
        produces<std::vector<reco::PFCandidate> >();
        muontype_=pfmuons;}
    else{
        muontype_=pfmuons; //just to avoid warnings
        throw std::runtime_error("muonRochesterCorrector: input muontype string not recognized");
    }

    //configure muonrochcorr
    muonRochCorr_=ztop::muonRochCorr(seed_);
    muonRochCorr_.setIsMC(isMC_);


}

void
muonRochesterCorrector::produce(edm::Event& evt, const edm::EventSetup& setup)
{


    if(muontype_ == recomuons){

        edm::Handle<std::vector<reco::Muon> > muons;
        evt.getByLabel(inputmuons_, muons);
        std::auto_ptr<std::vector<reco::Muon> > outputmuons(new std::vector<reco::Muon>() );

        for(std::vector<reco::Muon>::const_iterator muon=muons->begin();muon<muons->end();++muon){
            TLorentzVector p4in(muon->px(),muon->py(),muon->pz(),muon->energy());
            muonRochCorr_.correctP4(p4in,muon->charge());
            reco::Muon tmpmu=*muon;

            tmpmu.setP4(math::XYZTLorentzVector(p4in.Px(),p4in.Py(),p4in.Pz(),p4in.E()));
            if(debug) comparePts(muon,&tmpmu);
            outputmuons->push_back(tmpmu);
        }

        std::sort(outputmuons->begin(), outputmuons->end(), GreaterByPt<reco::Muon>());
        evt.put(outputmuons);
        return;
        ///...
    }
    if(muontype_ == pfmuons){
        edm::Handle<std::vector<reco::PFCandidate> > muons;
        evt.getByLabel(inputmuons_, muons);
        std::auto_ptr<std::vector<reco::PFCandidate> > outputmuons(new std::vector<reco::PFCandidate>() );

        for(std::vector<reco::PFCandidate>::const_iterator muon=muons->begin();muon<muons->end();++muon){
            TLorentzVector p4in(muon->px(),muon->py(),muon->pz(),muon->energy());
            muonRochCorr_.correctP4(p4in,muon->charge());
            reco::PFCandidate tmpmu=*muon;

            tmpmu.setP4(math::XYZTLorentzVector(p4in.Px(),p4in.Py(),p4in.Pz(),p4in.E()));
            if(debug) comparePts(muon,&tmpmu);
            outputmuons->push_back(tmpmu);
        }


        std::sort(outputmuons->begin(), outputmuons->end(), GreaterByPt<reco::PFCandidate>());
        evt.put(outputmuons);
        return;
        ///...
    }
    if(muontype_ == patmuons){

        //  std::vector<pat::Muon> * outmuons = new std::vector<pat::Muon>();
        edm::Handle<std::vector<pat::Muon> > muons;
        evt.getByLabel(inputmuons_, muons);
        std::auto_ptr<std::vector<pat::Muon> > outputmuons(new std::vector<pat::Muon>() );
        for(std::vector<pat::Muon>::const_iterator muon=muons->begin();muon<muons->end();++muon){
            TLorentzVector p4in(muon->px(),muon->py(),muon->pz(),muon->energy());
            muonRochCorr_.correctP4(p4in,muon->charge());
            pat::Muon tmpmu=*muon;

            tmpmu.setP4(math::XYZTLorentzVector(p4in.Px(),p4in.Py(),p4in.Pz(),p4in.E()));
            if(debug) comparePts(muon,&tmpmu);
            outputmuons->push_back(tmpmu);
        }


        std::sort(outputmuons->begin(), outputmuons->end(), GreaterByPt<pat::Muon>());
        evt.put(outputmuons);
        return;
        ///...
    }


} ////////////////////////////////////////////////////////////////////////

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE( muonRochesterCorrector );
