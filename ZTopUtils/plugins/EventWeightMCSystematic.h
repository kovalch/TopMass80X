///////////////////////////////////////////////////////////////////////////////
//
// 
// 15/05/06                   Alexander Grohsjean <alexander.grohsjean@desy.de>
// 
//
//                     print out to connect weight ID and systematic variation
//                                          filling of systematic event weight
//
//
///////////////////////////////////////////////////////////////////////////////
#ifndef EventWeightMCSystematic_h
#define EventWeightMCSystematic_h

#include <memory>
#include <string>
#include <iostream>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TopAnalysis/ZTopUtils/interface/version.h"

class EventWeightMCSystematic : public edm::EDProducer {
    
public:
    explicit EventWeightMCSystematic(const edm::ParameterSet&);
    ~EventWeightMCSystematic();
    
private:
    virtual void produce(edm::Event&, const edm::EventSetup&);
    #ifndef CMSSW_LEQ_5
    virtual void endRun(edm::Run const&, edm::EventSetup const&) override;  
    #endif 
    // ----------member data ---------------------------                                                 |
    edm::InputTag genEventInfoTag_;
    edm::InputTag lheEventInfoTag_;
    std::string weightID_;
    bool printLHE_;
};

#endif
