#ifndef muonRochesterCorrector_h
#define muonRochesterCorrector_h

#include <memory>
#include <string>
#include <iostream>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "../interface/muonRochCorr.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include <iostream>
/**
 * This small module can combine bools from the event content
 * It returns true if all (if any) sources from "mustBeTrue" are true
 * and all sources (if any) from "mustBeFalse" are false
 * Else it returns false
 *
 * This is NOT a filter itself so it does not stop events
 * from being processed
 */
class muonRochesterCorrector : public edm::EDProducer {

public:
    explicit muonRochesterCorrector(const edm::ParameterSet&);
    ~muonRochesterCorrector(){};

private:
    virtual void produce(edm::Event&, const edm::EventSetup&);
    template<class T,class U>
    void comparePts(T  a, U*b){
        if(debug){
            std::cout << a->pt() << "     " << b->pt();
            std::cout << std::endl;}
    }

private:
    enum muontypes{recomuons,pfmuons,patmuons};

    bool debug;
    edm::InputTag inputmuons_ ;
    bool isMC_;
    std::string muontypestring_;
    muontypes muontype_;
    int seed_;
    ztop::muonRochCorr  muonRochCorr_;



};



#endif
