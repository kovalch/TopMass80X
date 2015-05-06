///////////////////////////////////////////////////////////////////////////////
//
// 
// 15/05/06                   Alexander Grohsjean <alexander.grohsjean@desy.de>
// 
//
// currently assuming that meaning of weight IDs is known to user
// not enough info in ntuple to read it of
// if weight ID is not found, default MC weight is returned
// code requires LHEEventProduct("externalLHEProducer") and 
// GenEventInfoProduct("generator") 
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


class EventWeightMCSystematic : public edm::EDProducer {
    
public:
    explicit EventWeightMCSystematic(const edm::ParameterSet&);
    ~EventWeightMCSystematic();
    //    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions); 
    
private:
    virtual void produce(edm::Event&, const edm::EventSetup&);
    
    // ----------member data ---------------------------                                                 |
    const edm::ParameterSet parameterSet_;

};

#endif
