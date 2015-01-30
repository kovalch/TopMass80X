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
    double weight_ ;  

};

#endif
