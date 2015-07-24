#ifndef MyGenJetsProducer_h
#define MyGenJetsProducer_h

#include <memory>
#include <string>
#include <iostream>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"

class MyGenJetsProducer : public edm::EDProducer {

 public:
  explicit MyGenJetsProducer(const edm::ParameterSet&);
  ~MyGenJetsProducer();
  
 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

 private:
  edm::InputTag jetsSrc_, leptonSrc_;

};

#endif
