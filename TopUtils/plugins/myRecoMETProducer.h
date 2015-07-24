#ifndef MyRecoMETProducer_h
#define MyRecoMETProducer_h

#include <memory>
#include <string>
#include <iostream>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"

class MyRecoMETProducer : public edm::EDProducer {

 public:
  explicit MyRecoMETProducer(const edm::ParameterSet&);
  ~MyRecoMETProducer();
  
 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

 private:
  edm::InputTag mETSrc_;

};

#endif
