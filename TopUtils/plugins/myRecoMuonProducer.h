#ifndef MyRecoMuonProducer_h
#define MyRecoMuonProducer_h

#include <memory>
#include <string>
#include <iostream>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"

class MyRecoMuonProducer : public edm::EDProducer {

 public:
  explicit MyRecoMuonProducer(const edm::ParameterSet&);
  ~MyRecoMuonProducer();
  
 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

 private:
  edm::InputTag muonSrc_;

};

#endif
