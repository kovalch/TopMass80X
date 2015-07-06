#ifndef EventWeightPDF_h
#define EventWeightPDF_h

#include <memory>
#include <string>
#include <iostream>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


class EventWeightPDF : public edm::EDProducer {
    
public:
    explicit EventWeightPDF(const edm::ParameterSet&);
    ~EventWeightPDF();
    
private:
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endRun(edm::Run const&, edm::EventSetup const&) override;    
    // ----------member data ---------------------------                                                 |
    edm::InputTag lheEventInfoTag_;
    std::vector<std::string> PDFName_;
    std::vector<std::string> shortPDFName_;
    std::vector<std::string> nominalWeightID_;
    std::vector<std::string> beginWeightID_;
    std::vector<std::string> endWeightID_;
    bool printLHE_;
};

#endif
