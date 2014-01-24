#ifndef Boolcombiner_h
#define Boolcombiner_h

#include <memory>
#include <string>
#include <iostream>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

/**
 * This small module can combine bools from the event content
 * It returns true if all (if any) sources from "mustBeTrue" are true
 * and all sources (if any) from "mustBeFalse" are false
 * Else it returns false
 *
 * This is NOT a filter itself so it does not stop events
 * from being processed
 */
class Boolcombiner : public edm::EDProducer {

 public:
  explicit Boolcombiner(const edm::ParameterSet&);
  ~Boolcombiner(){};

 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

 private:
  bool debug;
  std::vector<edm::InputTag> mustbetrue_,mustbefalse_ ;
};



#endif
