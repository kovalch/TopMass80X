#ifndef MyGenElectronToMuonConverter_h
#define MyGenElectronToMuonConverter_h

#include <memory>
#include <string>
#include <iostream>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"

class MyGenElectronToMuonConverter : public edm::EDProducer {

 public:
  explicit MyGenElectronToMuonConverter(const edm::ParameterSet&);
  ~MyGenElectronToMuonConverter();
  
 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

 private:
  edm::InputTag electronSrc_;

};

#endif
